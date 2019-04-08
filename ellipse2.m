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
major = major / 2; % radii rather than diameters
minor = minor / 2;
clear sz pos

%% Set-Up
bw = zeros(m,n);

ind = 1:(m*n); % linear indices of all pixels in bw
[I,J] = ind2sub([m,n], ind);

boxmask = (abs(I-y0)<major) & (abs(J-x0)<major); % a bounding box enclosing the ellipse (should make the function slightly faster by only considering these pixels)

Ibox = I(boxmask);
Jbox = J(boxmask);

%% Transform coordinates:
Iellipse = Ibox - y0; % translating coordinates
Jellipse = Jbox - x0;

rotmat = [cosd(-phi),sind(-phi);-sind(-phi),cosd(-phi)]; % rotation matrix
IJellipse = [Iellipse',Jellipse'] * rotmat; % coordinates of box pixels in the ellipse's coordinate frame

%% Ellipse equation:
ellipsemask = (IJellipse(:,2).^2 ./ (major.^2)) + (IJellipse(:,1).^2 ./ (minor.^2)) < 1; % binary mask for Ibox and Jbox, indicating which pixels lie inside the ellipse

%% Draw ellipse:
Iellipse = Ibox(ellipsemask');
Jellipse = Jbox(ellipsemask');
ellilpseIdx = sub2ind([m,n], Iellipse, Jellipse); % linear indices of ellipse pixels in bw
bw(ellilpseIdx) = 1;

end