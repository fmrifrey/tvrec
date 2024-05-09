function P = L_adj(x)
% L_adj() computes the adjoint operation P = L'x as described in section 4
%   of Beck, Amir, and Marc Teboulle. “Fast Gradient-Based Algorithms for
%   Constrained Total Variation Image Denoising and Deblurring Problems.”
%
% written by David Frey (djfrey@umich.edu) and Tao Hong (tahong@umich.edu)
%
% inputs:
%     x             matrix of any size
%                       
% outputs:
%     P             matrix set {p,q,...} of difference matrices along each
%                       dimension; size() = {ndims(x),1}
%

    % get size
    sz = size(x);
    nd = ndims(x);
    if (nd==2 && sz(2)==1)
        nd = 1;
    end

    % initialize cell array of finite diff matrices
    P = cell(nd,1);
    
    % loop through dimensions
    for d = 1:nd
        P{d} = -diff(x,1,d);
    end

end

