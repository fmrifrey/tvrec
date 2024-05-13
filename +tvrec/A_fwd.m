function b = A_fwd(x,F,ncoils,parallelize)
% A_adj() computes the forward operation: b = A(x) for NUFFT-SENSE
%   reconstruction
%
% written by David Frey (djfrey@umich.edu) and Tao Hong (tahong@umich.edu)
%
% inputs:
%     x             image space data; size() = [N(:)', Nt]
%     F             cell array of system matrices for each time point
%                       (should not include density compensation);
%                       size() = {Nt, 1}
%     ncoils        number of coils in output data
%     paralellize   option to run frame-wise reconstructions in parallel
%                       (0 or 1)
% 
% outputs:
%     b             kspace data; size() = [Nk, Nt, Nc]
%
    
    % set default parallelize
    if nargin < 4 || isempty(parallelize)
        parallelize = 0;
    end
    
    % get sizes
    Nk = F{1}.odim(1);
    Nt = length(F);
    Nc = ncoils;

    % initialize b
    b = zeros(Nk,Nt,Nc);
    x = reshape(x,prod(F{1}.idim),Nt);
    parfor (n = 1:Nt, 100*parallelize) % loop through time points
        A = F{n};
        b(:,n,:) = reshape(A*reshape(x(:,n),F{1}.idim),[],1,ncoils); % calculate NUFFT of image at time point n
    end

end