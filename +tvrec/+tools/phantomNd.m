function [P,E] = phantomNd(N,E)
% [P,E] = phantomNd(N,E) computes a 2d-4d ellipsoid-based digital phantom
%   image, P given image dimensions N and phantom specs, E
%
% N is an [Nd x 1] vector containing grid sizes of image P. 4th
%     dimension is assumed to be time (2d timeseries data should be passed
%     as [Nx x Ny x 1 x Nt])
%
% E is a matrix representing a user-defined phantom, where each row
%     of the matrix E specifies an ellipsoid in the image. E is of
%     size [Nellipsoids x 10 x Nt], with each column containing a
%     different parameter for the ellipsoids:
%         Column 1:  A      the additive intensity value of the ellipsoid
%         Column 2:  a      the length of the x semi-axis of the ellipsoid 
%         Column 3:  b      the length of the y semi-axis of the ellipsoid
%         Column 4:  c      the length of the z semi-axis of the ellipsoid
%           (value not used in 2d case)
%         Column 5:  x0     the x-coordinate of the center of the ellipsoid
%         Column 6:  y0     the y-coordinate of the center of the ellipsoid
%         Column 7:  z0     the z-coordinate of the center of the ellipsoid
%           (value not used in 2d case)
%         Column 8:  phi    phi Euler angle (in degrees) (z rotation 2)
%         Column 9:  theta  theta Euler angle (in degrees) (x rotation)
%         Column 10: psi    psi Euler angle (in degrees) (z rotation 1)
%           (value not used in 2d case)
%
% E can also be passed as a string that specifies the type of default head
%     phantom to generate. Valid values are:
%      'Shepp-Logan'            A test image used widely by researchers in
%                               tomography
%      'Modified Shepp-Logan'   (default) A variant of the Shepp-Logan
%                               phantom in which the contrast is improved
%                               for better visual perception.
%       'Yu-Ye-Wang'            A phantom described in Yu H, Ye Y, Wang G,
%                               Katsevich-Type Algorithms for Variable
%                               Radius Spiral Cone-Beam CT
%
% If Nt > 1, but size(E,3) == 1, random motion of the ellipsoids will be
%     simulated
%
% Credit to Matthias Christian Schabel (matthias@stanfordalumni.org)
%     for original phantom3d.m function
%

% set default phantom
if nargin < 2 || isempty(E)
    E = 'Modified Shepp-Logan';
end

% load default head phantom
if ischar(E)
    switch lower(E)
        case 'shepp-logan'
            E = shepp_logan;
        case 'modified shepp-logan'
            E = modified_shepp_logan;
        case 'yu-ye-wang'
            E = yu_ye_wang;
        otherwise
            error('invalid default phantom: %s',E);
    end
end
    
% simulate random noise along time
N = padarray(N(:),[4-length(N),0],1,'post');
if N(4) > 1 && size(E,3) == 1
    E = repmat(E,[1,1,1,N(4)]);
    E(:,5:end,2:end) = E(:,5:end,2:end) + ...
        ... % vary a-c, x0-z0, and phi-psi by 5% over time
        0.1*E(:,5:end,2:end).*(2*rand(size(E,1),6,N(4)-1)-1);
elseif N(4) ~= size(E,3)
    error('size(E,3) must be == Nt, or 1 to simulate random noise');
end

% make image grid points
x = cell(3,1);
for i = 1:3
    if N(i) > 1
        x{i} = linspace(-1,1,N(i));
    else
        x{i} = zeros(1,N(i));
    end
end
[X,Y,Z] = ndgrid(x{:});
r = [X(:),Y(:),Z(:)]';

% initialize vectorized image array
P = zeros(size(r,2),N(4));

for nt = 1:N(4) % loop through time points
    for ne = 1:size(E,1) % loop through ellipsoids

        % get properties of current ellipse
        rho = E(ne,1,nt); % intensity
        D = diag(E(ne,2:4,nt)); % stretch
        rc = E(ne,5:7,nt)'; % center
        R = eul2rotm(pi/180*E(ne,10:-1:8,nt),'ZYZ'); % rotation
        
        % determine ellipsoid ROI and add amplitude
        ROI = vecnorm(D^-1*R'*(r - rc), 2, 1) <= 1;
        P(ROI(:),nt) = P(ROI(:),nt) + rho;

    end
end

% reshape the image
P = reshape(P,N');

end

function e = shepp_logan

e = modified_shepp_logan;
e(:,1) = [1 -.98 -.02 -.02 .01 .01 .01 .01 .01 .01];

end
      
function e = modified_shepp_logan
%
%   This head phantom is the same as the Shepp-Logan except 
%   the intensities are changed to yield higher contrast in
%   the image.  Taken from Toft, 199-200.
%      
%         A      a     b     c     x0      y0      z0    phi  theta    psi
%        -----------------------------------------------------------------
e =    [  1  .6900  .920  .810      0       0       0      0      0      0
        -.8  .6624  .874  .780      0  -.0184       0      0      0      0
        -.2  .1100  .310  .220    .22       0       0    -18      0     10
        -.2  .1600  .410  .280   -.22       0       0     18      0     10
         .1  .2100  .250  .410      0     .35    -.15      0      0      0
         .1  .0460  .046  .050      0      .1     .25      0      0      0
         .1  .0460  .046  .050      0     -.1     .25      0      0      0
         .1  .0460  .023  .050   -.08   -.605       0      0      0      0
         .1  .0230  .023  .020      0   -.606       0      0      0      0
         .1  .0230  .046  .020    .06   -.605       0      0      0      0 ];
       
end       

function e = yu_ye_wang
%
%   Yu H, Ye Y, Wang G, Katsevich-Type Algorithms for Variable Radius Spiral Cone-Beam CT
%      
%         A      a     b     c     x0      y0      z0    phi  theta    psi
%        -----------------------------------------------------------------
e =    [  1  .6900  .920  .900      0       0       0      0      0      0
        -.8  .6624  .874  .880      0       0       0      0      0      0
        -.2  .4100  .160  .210   -.22       0    -.25    108      0      0
        -.2  .3100  .110  .220    .22       0    -.25     72      0      0
         .2  .2100  .250  .500      0     .35    -.25      0      0      0
         .2  .0460  .046  .046      0      .1    -.25      0      0      0
         .1  .0460  .023  .020   -.08    -.65    -.25      0      0      0
         .1  .0460  .023  .020    .06    -.65    -.25     90      0      0
         .2  .0560  .040  .100    .06   -.105    .625     90      0      0
        -.2  .0560  .056  .100      0    .100    .625      0      0      0 ];
       
end