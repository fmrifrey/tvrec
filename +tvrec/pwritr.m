function L = pwritr(A,At,sz,tol)
% pwritr() estimates the Lipschitz constant, L for a given system operator
%   A using the power iteration method
%
% written by David Frey (djfrey@umich.edu) and Tao Hong (tahong@umich.edu)
%
% inputs:
%      A            forward system operator (function handle)
%      At           adjoint system operator (function handle)
%      sz           size of the input image; size() = [nd,1]
%      tol          convergence tolerance (percent error)
%
% outputs:
%      L            estimated Lipschitz constant
%

    % set default tolerance
    if nargin < 4 || isempty(tol)
        tol = 1; % 1% error
    end

    % initialize
    b_k = randn(sz);
    Ab_k = At(A(b_k));
    b_k_1 = b_k;
    norm_b_k = norm(Ab_k(:),'fro');
    
    while 1
    
        % calculate new norm
        b_k = Ab_k/norm_b_k;
        Ab_k = (At(A(b_k)));
        norm_b_k_1 = norm(Ab_k(:),'fro');
        
        % check for convergence
        if 100*norm(b_k_1(:)-b_k(:),'fro')/norm(b_k(:),'fro') <= tol
            break
        else
            norm_b_k = norm_b_k_1;
            b_k_1 = b_k;
        end    
    
    end
    b = b_k;
    
    % calculate Lipschitz constant
    L = real(b(:)'*Ab_k(:)/(b(:)'*b(:)));
    
end

