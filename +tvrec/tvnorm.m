function tv = tvnorm(x,type)
% tvnorm() calculates the discrete total variation semi-norm, described by
%   the equations following equation 2.2 in Beck, Amir, and Marc Teboulle.
%   “Fast Gradient-Based Algorithms for Constrained Total Variation Image
%   Denoising and Deblurring Problems.”
%
% written by David Frey (djfrey@umich.edu) and Tao Hong (tahong@umich.edu)
%
% input
%     x     image matrix; size() = [N(:)', Nt]
%     type  type of norm ('l1' or 'iso')
% outputs:
%     tv    value of discrete total variation semi-norm
%

    % set default type
    if nargin < 2 || isempty(type)
        type = 'l1';
    end

    switch type
        case 'iso'
            tv = tvnorm_iso(x);
        case 'l1'
            tv = tvnorm_l1(x);
        otherwise
            error('invalid type: %s',type);
    end

end

function tv = tvnorm_iso(x)

    % get size
    sz = size(x);
    nd = ndims(x);
    if (nd==2 && sz(2)==1)
        nd = 1;
    end

    % calculate projection operators
    P = tvrec.L_adj(x);

    % loop through dimensions
    D = zeros(sz); % sum of squares matrix
    for d = 1:nd
        % add its square to D
        D = D + P{d}.^2;
    end

    % calculate norm
    tv = sum(sqrt(D),'all');

end

function tv = tvnorm_l1(x)

    % get size
    sz = size(x);
    nd = ndims(x);
    if (nd==2 && sz(2)==1)
        nd = 1;
    end

    % calculate projection operators
    P = tvrec.L_adj(x);

    % loop through dimensions
    tv = 0;
    for d = 1:nd
        % add l1 norm
        tv = tv + sum(abs(P{d}),'all');
    end
    
end