function L = pwritr(A,At,sz,tol,show)
% pwritr() estimates the Lipschitz constant, L for a given system operator
%   A using the power iteration method
%
% written by David Frey (djfrey@umich.edu) and Tao Hong (tahong@umich.edu)
%
% inputs:
%      A            forward system operator (function handle)
%      At           adjoint system operator (function handle)
%      sz           size of the input image; size() = [nd,1]
%      tol          convergence tolerance
%      show         option to show power iteration convergence (0 or 1)
%
% outputs:
%      L            estimated Lipschitz constant
%

    % set default tolerance
    if nargin < 4 || isempty(tol)
        tol = 1e-2;
    end
    
    % set default tolerance
    if nargin < 5 || isempty(show)
        show = 1;
    end
    
    % initialize
    b_k = randn(sz);
    Ab_k = At(A(b_k));
    b_k_1 = b_k;
    norm_b_k = norm(Ab_k(:),'fro');
    y = [];
    
    while 1
    
        % calculate new norm
        b_k = Ab_k/norm_b_k;
        Ab_k = (At(A(b_k)));
        norm_b_k_1 = norm(Ab_k(:),'fro');
        
        % check for convergence
        if norm(b_k_1(:)-b_k(:),'fro')/norm(b_k(:),'fro') <= tol
            break
        end    
        
        if show
            y = [y; norm(b_k_1(:)-b_k(:),'fro')/norm(b_k(:),'fro')];
            cfigopen('power iteration');
            plot(y), hold on
            yline(tol,'--r','Label',sprintf('tol = %.1f%%',tol)); hold off
            xlabel('# iterations');
            ylabel('norm error');
        end
        
        norm_b_k = norm_b_k_1;
        b_k_1 = b_k;
    
    end
    b = b_k;
    
    % calculate Lipschitz constant
    L = real(b(:)'*Ab_k(:)/(b(:)'*b(:)));
    
end

