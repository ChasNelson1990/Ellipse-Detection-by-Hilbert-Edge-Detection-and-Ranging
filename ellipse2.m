function bw = ellipse2(sz,pos,major,minor,phi)
%% ellipse2, a function to produce a binary image containing a single ellipse
%
% ellipse2(sz,pos,major,minor,phi)
%
% Details   This function produces a binary image containing a single
%           filled ellipse with the given parameters.
% Inputs    sz - desired size of image, if scalar image is square, if a
%           vector image is sz(1) by sz(2)
%           pos - centroid position within the image, if scalar ellipse is
%           centere at [pos,pos], if a vector at [pos(1),pos(2)]
%           major - major axis length (full) of the ellipse in pixels
%           minor - minor axis length (full) of the ellipse in pixels
%           phi - rotation of ellipse, i.e. angle going clockwise from the
%           x axis of the image to the major axis of the ellipse
% Outputs   bw - binary image
%
% Examples:
% bw = ellipse2(20,10,6,2), creates an image of size 20x20 with an ellipse
% centred at x=10, y=10 of size 6 (major) by 2 (minor) and orientated along
% the x-axis
% bw = ellipse2([100,50],[20,40],10,5,57), creates an image of size 100x50
% with an ellipse centred at x=20, y=40 of size 10 (major) by 5 (minor)
% nd orientated along 57 degrees to the x-axis
% the x-axis
%
% Copyright 2016 Carl J. Nelson, Durham University, UK
%
% License   See included <a href="./LICENSE/">file</a> or visit
%           <a href="https://github.com/ChasNelson1990/...
%           	Ellipse-Detection-by-Hilbert-Edge-Detection-and-Ranging/">
%			The GitHub Repository</a>

%% Inputs
if nargin<5
    phi = 0;
end
if length(sz)==2
    m = sz(1); n = sz(2);
else
    m = sz; n = sz;
end
x0 = pos(1); y0 = pos(2);
clear sz pos

%% Set-Up
bw = zeros(m,n);

%% Circumference & Parameterisation in Theta
h = (major-minor)^2/(major+minor)^2;
circ = pi*(major+minor)*(1+((3*h)/(10+sqrt(4-(3*h)))));%Approximate (Ramanujan, 1914)
theta = linspace(0,2*pi,100*round(circ));

%% Calculate Ellipse Outline Pixels
x = x0 + 0.5*(major*cos(theta)*cosd(-phi) - minor*sin(theta)*sind(-phi));
y = y0 + 0.5*(major*cos(theta)*sind(-phi) + minor*sin(theta)*cosd(-phi));

%% Boundary Issues
x(x<1) = 1; y(y<1) = 1;
x(x>n) = n; y(y>m) = m;

%% Draw Outline
idx = sub2ind(size(bw),round(y),round(x));
bw(idx) = 1;

%% Fill Ellipse
bw = imfill(bw,'holes');

end
