function isadj = testadj(A,At,xtmp,niter,tol)
% testadj() calculates the ratio <A(x),y>/<A'(y),x> for niter random
%   inputs, x to test for adjointness of the forward and adjoint
%   operators, A and A'
%
% written by David Frey (djfrey@umich.edu)
%
% inputs:
%     A             forward system operator (function handle)
%     At            apparent adjoint system operator (function handle)
%     xtmp          template input matrix of same size as x (can be a cell
%                       array of many matrices if consistent with operator)
%     niter         number of iterations to perform
%     tol           mean percent error tolerance to return iadj = 1
% 
% outputs:
%     isadj         binary claiming if A and A' are adjoint (0 or 1); if no
%                       output is returned, function will show a box plot
%                       of the adjoint ratios
%

if nargin<3 || isempty(xtmp)
    xtmp = zeros(64);
end

if nargin<4 || isempty(niter)
    niter = 500;
end

if nargin < 5 || isempty(tol)
    tol = 1e-2;
end

ratios = zeros(niter,1);
for itr = 1:niter

    % generate random x and compute Ax
    if iscell(xtmp) % for cases that x is a cell
        x = cell(size(xtmp));
        for i = 1:length(x)
            x{i} = randn(size(xtmp{i}));
        end
    else
        x = randn(size(xtmp));
    end
    Ax = A(x);
    
    % generate random y and compute A'y
    if iscell(Ax) % for cases that y is a cell
        y = cell(size(Ax));
        for i = 1:length(y)
            y{i} = randn(size(Ax{i}));
        end
    else
        y = randn(size(Ax));
    end
    Aty = At(y);
    
    % calculate numerator <Ax,y>
    if iscell(y) % for cases that y is a cell
        numer = 0;
        for i = 1:length(y)
            numer = numer + Ax{i}(:)'*y{i}(:);
        end
    else
        numer = Ax(:)'*y(:);
    end
    
    % calculate denominator <x,A'y>
    if iscell(x) % for cases that y is a cell
        denom = 0;
        for i = 1:length(x)
            denom = denom + x{i}(:)'*Aty{i}(:);
        end
    else
        denom = x(:)'*Aty(:);
    end
    
    % calculate <Ax,y>/<x,A'y>
    ratios(itr) = abs(numer/denom);

end

% plot histogram of ratios
if nargout < 1
    tvrec.tools.cfigopen('adjointness test')
    boxplot(ratios)
    ylabel("|<Ax,y>/<x,A'y>|")
    xticks([])
else
    isadj = 100 * abs(mean(ratios) - 1) < tol;
end

end

