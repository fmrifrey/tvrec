function b = A_fwd(x,F,smap,parallelize)
% x = image {[N x N] x Nt}
% F = NUFFT operators {[Nk x N^2] x Nt}
% smap = sensitivity map [N x N]
% parallelize = option to parallelize the frame loop (0 or 1)
    
    % set default parallelize
    if nargin < 4 || isempty(parallelize)
        parallelize = 0;
    end
    
    % get sizes
    Nk = F{1}.odim(1);
    Nt = length(F);
    Nc = size(smap,length(F{1}.idim)+1);

    % initialize b
    b = zeros(Nk,Nc,Nt);
    x = reshape(x,prod(F{1}.idim),Nt);
    parfor (n = 1:Nt, 100*parallelize) % loop through time points
        A = F{n};
        if ~isempty(smap)
            A = Asense(A,smap); % sense
        end
        b(:,:,n) = A*reshape(x(:,n),F{1}.idim); % calculate NUFFT of image at time point n
    end

end