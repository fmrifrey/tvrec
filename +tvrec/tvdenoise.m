function [Y,P] = tvdenoise(v,lam,P,type,niter)
% v = noisey image
% lam = lagrange multiplier
% P = cell array of projection operators {Nd x 1}
% type = type of tv semi-norm ('iso' or 'l1')
% niter = max number of iterations

    % get size
    sz = size(v);
    nd = ndims(v);
    if (nd==2 && sz(2)==1)
        nd = 1;
    end

    % set default P
    if isempty(P)
        for d = 1:nd
            szp = sz;
            szp(d) = szp(d)-1;
            P{d} = eps()*ones(szp);
        end
    end

    % initialize
    R = P;
    D = zeros(sz);
    t_k_1 = 1;
    ep = 1e-5;
    
    % loop through iterations
    for i = 1:niter
    
        % store old values
        D_old = D;
        P_old = P;
        t_k = t_k_1;

        % compute gradient of objective fun
        D = v - lam * tvrec.L_fwd(R);
       
        Q = tvrec.L_adj(D);

        % take a step towards negative of the gradient
        for d = 1:nd
            P{d} = R{d} + 1/(4*nd*lam)*Q{d};
        end
        
        % calculate projection
        switch type
            case 'iso'
                P = proj_iso(P,sz,nd);
            case 'l1'
                P = proj_l1(P,nd);
            otherwise
                error('invalid type: %s',type);
        end
        
        % calculate t_{k+1}
        t_k_1 = (1 + sqrt(1+4*t_k^2))/2;

        % update R
        for d = 1:nd
            R{d} = P{d} + (t_k-1)/t_k_1 * (P{d} - P_old{d});
        end

        % calculate residual
        res = norm(D(:) - D_old(:),'fro') / norm(D(:),'fro');
        if (res < ep)
            break
        end
        
    end

    % calculate Y
    Y = v - lam * tvrec.L_fwd(P);

end

function P = proj_iso(P,sz,nd)

    % loop through dimensions
    A = zeros(sz); % root sum of squares
    for d = 1:nd
        padsz = zeros(nd,1);
        padsz(d) = 1;
        % add square of P to A
        A = A + padarray(P{d},padsz,0,'post').^2;
    end

    % take root of sum of squares
    A = sqrt(max(A,1));

    % loop through dimensions
    for d = 1:nd
        % divide P by rsos
        P{d} = P{d}./A;
    end

end

function P = proj_l1(P,nd)
    
    % loop through dimensions
    for d = 1:nd
        % divide P by absolute max
        Pd = P{d};
        P{d} = P{d}./max(abs(Pd),1);
    end

end