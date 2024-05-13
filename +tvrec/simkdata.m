function [kdata,smap,x_gt] = simkdata(klocs,N,fov,varargin)
% simkdata() simulates kspace data at sampling locations for digital
%   ellipsoidal phantom data
%
% written by David Frey (djfrey@umich.edu) and Tao Hong (tahong@umich.edu)
%
% inputs:
%     klocs         kspace sampling locations; size() = [Nk, Nt, Nd]
%     N             image dimensions; size() = [Nd, 1]
%     fov           field of view (in inverse units of klocs);
%                       size() = [Nd, 1]
%     'ncoils'      number of coils to simulate
%     'snr'         signal to noise ratio for additive white gaussian
%                       noise (dB)
%     'show'        option to show the sampling (2d data only); (0 or 1)
% 
% outputs:
%      kdata        kspace data at sampling locs; size() = [Nk, Nt, Nc]
%      smap         simulated sensitivity maps; size() = [N(:)', Nc]
%      x_gt         ground truth digital phantom image;
%                       size() = [N(:)', Nt]
%

    % define defaults
    defaults = struct( ...
        'ncoils', 4, ... % number of coils to simulate
        'snr', inf, ... % signal to noise ratio (see agwn.m)
        'show', 1 ... % show the ground truth image and sampling scheme
        );

    % parse arguments
    arg = vararg_pair(defaults,varargin);

    N = N(:)';
    fov = fov(:)';
    nd = size(N,2);
    Nt = size(klocs,2);
    
    % generate ground truth phantom image & sensitivity maps
    x_gt = tvrec.tools.phantomNd([N,Nt]);
    switch nd
        case 2
            smap = ir_mri_sensemap_sim('nx', N(1), 'ny', N(2),  ...
                'dx', fov(1)/N(1), 'dy', fov(2)/N(2), ...
                'rcoil', 10, 'ncoil', arg.ncoils, 'chat', 0);
        case 3
            smap = ir_mri_sensemap_sim('nx', N(1), 'ny', N(2), 'nz', N(3), ...
                'dx', fov(1)/N(1), 'dy', fov(2)/N(2), 'dz', fov(3)/N(3), ...
                'rcoil', 10, 'ncoil', arg.ncoils, 'chat', 0);
        otherwise
            error('nd must be either 2 or 3 (2d or 3d)');
    end
    if arg.ncoils == 1
        smap = ones(size(smap));
    end
    
    % create system matrices
    A = cell(Nt,1);
    for i = 1:Nt
        % reformat for NUFFT
        omega = 2*pi*fov./N.*squeeze(klocs(:,i,1:nd));
        
        % generate the fwd operator using NUFFT
        nufft_args = {N, 6*ones(1,nd), 2*N, N/2, ...
            'table', 2^10, 'minmax:kb'};
        A{i} = Gnufft(true(N),cat(2,{omega},nufft_args)); % NUFFT
        if i==1
            A1 = A{i};
        end
        if arg.ncoils > 1
            A{i} = Asense(A{i},smap);
        end
    end
    
    % NUFFT the object - inverse crime k-space signal estimation
    kdata = tvrec.A_fwd(x_gt,A,0);
    kdata = awgn(kdata,arg.snr); % add white gaussin noise
    
    % create a figure showing k-space sampling pattern
    if arg.show && nd == 2
        tvrec.tools.cfigopen('simpdata(): k-space sampling')
        
        subplot(1,3,1)
        xx = fov(1)*linspace(-0.5,0.5,N(1));
        yy = fov(2)*linspace(-0.5,0.5,N(2));
        im(xx,yy,abs(x_gt));
        xlabel('x (cm)');
        ylabel('y (cm)');
        title('ground truth')
        colormap parula
        
        subplot(1,3,2)
        for i = 1:size(klocs,2)
            kdata2 = reshape(kdata(:,:,1),size(klocs,1),size(klocs,2),[]);
            col = mean(abs(kdata2(:,i,:)),3)./max(abs(kdata2(:,i,:)),[],'all');
            x = klocs(:,i,1)';
            y = klocs(:,i,2)';
            z = zeros(1,size(klocs,1));
            c = col';
            surf([x;x],[y;y],[z;z],[c;c], ...
                'facecol', 'no', ...
                'edgecol', 'interp', ...
                'linew', 2 ...
                );
            hold on
        end
        hold off
        view(2)
        set(gca,'color','k');
        xlabel('kx (cm^{-1})');
        ylabel('ky (cm^{-1})');
        title('k-space sampling');
        
        subplot(1,3,3)
        [~,tmpw] = tvrec.pmdcf(A1);
        im(xx,yy,abs(tvrec.A_adj(kdata,A(1),{tmpw},0)));
        title("A'y reconstruction");
        xlabel('x (cm)');
        ylabel('y (cm)');
        colormap parula
        
        drawnow
    elseif arg.show
        warning('figures for nd > 2 are not yet supported, sorry!\n');
    end

end

