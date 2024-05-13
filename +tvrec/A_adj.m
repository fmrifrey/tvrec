function x = A_adj(b,F,w,parallelize)
% A_adj() computes the adjoint operation: x = A'(b) for NUFFT-SENSE
%   reconstruction
%
% written by David Frey (djfrey@umich.edu) and Tao Hong (tahong@umich.edu)
%
% inputs:
%     b             kspace data; size() = [Nk, Nt, Nc]
%     F             cell array of system matrices for each time point
%                       (should not include density compensation or
%                       sensitvity encoding); size() = {Nt, 1}
%     w             cell array of sampling weighting arrays for each time
%                       point; size() = {Nt, 1}
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
    
    % get number of frames
    Nt = length(F);

    % initialize x
    x = zeros(prod(F{1}.idim),Nt);
    
    parfor (n = 1:Nt, 100*parallelize) % loop through time points
        A = F{n};
        xn = A'*squeeze(w{n}.*b(:,n,:)); % calculate inverse NUFFT of ksignal at time point n
        x(:,n) = reshape(xn,[],1);
    end
    x = reshape(x,[F{1}.idim,Nt]);

end