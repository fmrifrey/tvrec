function b = A_fwd(x,F,smap,parallelize)
% A_adj() computes the forward operation: b = A(x) for NUFFT-SENSE
%   reconstruction
%
% written by David Frey (djfrey@umich.edu) and Tao Hong (tahong@umich.edu)
%
% inputs:
%     b             kspace data; size() = [Nk, Nt, Nd]
%     F             cell array of system matrices for each time point
%                       (should not include density compensation or
%                       sensitvity encoding); size() = {Nt, 1}
%     smap          sensitivity maps, leave empty for no sensitivity
%                       encoding; size() = [N(:)', Nc]
%     paralellize   option to run frame-wise reconstructions in parallel
%                       (0 or 1)
% 
% outputs:
%      x            image space data; size() = [N(:)', Nt]
%
    
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