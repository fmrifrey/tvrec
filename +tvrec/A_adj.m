function x = A_adj(b,F,W,smap,parallelize)
% b = kspace data [Nk x Nc x Nt]
% F = NUFFT operators {[Nk x N^Nd] x Nt}
% W = density weighting matrix {[Nk x Nk] x Nt}
% smap = sensitivity map [N x N]
% parallelize = option to parallelize the frame loop (0 or 1)
    
    % set default parallelize
    if nargin < 4 || isempty(parallelize)
        parallelize = 0;
    end
    
    % get number of frames
    Nt = length(F);

    % initialize x
    x = zeros(prod(F{1}.idim),Nt);
    
    parfor (n = 1:Nt, 100*parallelize) % loop through time points
        A = F{n};
        A = W{n}*A; % density compensation
        if ~isempty(smap)
            A = Asense(A,smap); % sense
        end
        xn = A'*b(:,:,n); % calculate inverse NUFFT of ksignal at time point n
        x(:,n) = reshape(xn,[],1);
    end
    x = reshape(x,[F{1}.idim,Nt]);

end