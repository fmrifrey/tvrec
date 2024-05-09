function [W,w] = pmdcf(F,itrmax)
% pmdcf() computes the density compensation weighting matrix W for a given
%   NUFFT fatrix object using the methods described in Pipe, James G., and
%   Padmanabhan Menon. “Sampling Density Compensation in MRI: Rationale and
%   an Iterative Numerical Solution”.
%
% written by David Frey (djfrey@umich.edu)
%
% inputs:
%      F            NUFFT fatrix (Gnufft or Gmr object from MIRT)
%      itrmax       number of iterations to perform
% outputs:
%      W            density weighting fatrix
%      w            sample desnity weights
%

    % set default itrmax
    if nargin < 2 || isempty(itrmax)
        itrmax = 15;
    end
    
    % initialize weights to 1 (psf)
    w = ones(size(F,1),1);
    
    % loop through iterations
    for i = 1:itrmax
        
        % Pipe algorithm: W_{i+1} = W_{i} / (G * (G' * W_{i}))
        d = real( F.arg.st.interp_table(F.arg.st, ...
            F.arg.st.interp_table_adj(F.arg.st, w) ) );
        w = w ./ d;
        w = w / norm(w,1);
        W = Gdiag(w);
        
    end
    
end