function [x_star, cost, x_set, L] = tvrec_nufft(klocs,kdata,N,fov,varargin)
% tvrec_nufft() reconstructs image x given kspace sampling locations and
%   corresponding kspace data using total variation regularized iterative
%   SENSE reconstruction
%
% written by David Frey (djfrey@umich.edu) and Tao Hong (tahong@umich.edu)
%
% if no inputs are passed, test case (written at bottom of code) will run
%
% inputs:
%     klocs         kspace sampling locations; size() = [Nk, Nt, Nd]
%     kdata         kspace data at sampling locs; size() = [Nk, Nt, Nc]
%     N             image dimensions; size() = [Nd, 1]
%     fov           field of view (in inverse units of klocs);
%                       size() = [Nd, 1]
%     'lam'         lagrange multiplier for TV regularization;
%                       size() = [1,1]
%     'L'           Lipschitz constant (leave empty for power iteration
%                       estimation); size() = [1,1]
%     'type'        type of total variation cost ('l1' or 'iso')
%     'niter'       number of iterations; size() = [1,1]
%     'smap'        sensitivity maps, leave empty to coil compress;
%                       size() = [N(:)', Nc]
%     'paralellize' option to run frame-wise reconstructions in parallel
%                       (0 or 1)
%     'denscomp'    option to weigh kspace points by sampling density using
%                       pipe & menon dcf (0 or 1)
%     'show'        option to show the reconstruction as it iterates
%                       (0 or 1)
% 
% outputs:
%      x_star       final iteration estimate image; size() = [N(:)', Nt]
%      cost         vector containing cost at each iteration;
%                       size() = [niter, 1]
%      x_set        set of images at each iteration; size() = {niter,1};
%      L            Lipschitz constant (for recycling)
%

    % run test case if no inputs are passed
    if nargin < 1
        [x_star, cost, x_set] = test_tvrec_nufft;
        return
    end

    % define defaults
    defaults = struct( ...
        'lam', 0, ...
        'L', [], ...
        'type', 'l1', ...
        'niter', 100, ...
        'smap', [], ...
        'parallelize', 0, ...
        'denscomp', 1, ...
        'show', 0 ...
        );

    % parse arguments
    arg = vararg_pair(defaults,varargin);
    
    % simulate data if none is passed
    if isempty(kdata)
        warning('no data passed, creating phantom simulation...');
        [kdata,arg.smap] = tvrec.simkdata(klocs,N,fov,'show',arg.show);
    end

    % get sizes
    Nt = size(klocs,2);
    Nd = length(N);
    Nc = size(kdata,3);
    
    % convert N and fov to row vectors
    N = N(:)';
    fov = fov(:)';

    % check SENSE map
    if isempty(arg.smap) && Nc > 1
        warning('sense map is empty, compressing data to 1 coil...');
        kdata = ir_mri_coil_compress(kdata,'ncoil',1);
        Nc = 1;
    end
    
    % create a cell array of forward NUFFT operators for each time pt
    F = cell(Nt,1);
    w = cell(Nt,1);
    nufft_args = {N, 6*ones(Nd,1), 2*N, N/2, 'table', 2^10, 'minmax:kb'};
    
    % create nufft matrices for each frame
    parfor (n = 1:Nt, 100*arg.parallelize)
        omega = 2*pi*fov./N.*squeeze(klocs(:,n,1:Nd));
        F{n} = Gnufft(true(N),cat(2,{omega},nufft_args)); % NUFFT
        [~,w{n}] = tvrec.pmdcf(F{n}); % density compensation
        if Nc > 1 % sensitivity encoding
            F{n} = Asense(F{n},arg.smap);
        end
    end
    
    % define the forward and adjoint operators A and At
    A = @(x) tvrec.A_fwd(x,F,Nc,arg.parallelize);
    At = @(b) tvrec.A_adj(b,F,w,arg.parallelize);

    % calculate default L using power iteration
    if isempty(arg.L) && arg.niter > 0
        arg.L = tvrec.pwritr(A,At,[N(:)',Nt]);
    elseif isempty(arg.L)
        arg.L = 0;
    end

    % recon the data
    [x_star,cost,x_set] = tvrec.tvdeblur(A,At,kdata, ...
        'lam', arg.lam, ...
        'L', arg.L, ...
        'type', arg.type, ...
        'niter', arg.niter, ...
        'show', arg.show);
    
    % return L
    L = arg.L;

end

function [x_star, cost, x_set] = test_tvrec_nufft()
% test case for tvrec_nufft - radial undersampling
    
    % set parameters
    N = [256,256];
    fov = [24,24];
    nspokes = 32;

    % initialize radial sampling trajectory
    klocs = zeros(N(1),nspokes,2);
    spoke0 = [1;0;0]*N(1)/fov(1)/2.*linspace(-1,1,N(1));
    for i = 1:nspokes
        % transform initial spoke
        spokei = eul2rotm((i-1)*pi/nspokes*[1,0,0])*spoke0;
        klocs(:,i,:) = permute(spokei(1:2,:),[2,3,1]);
    end
    klocs = reshape(klocs,[],1,2); % combine spokes into 1 frame
    
    % run the reconstruction with simulated data
    [x_star, cost, x_set] = tvrec_nufft(klocs,[],N,fov,'lam',0.1,'show',1);

end