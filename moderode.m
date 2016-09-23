function B = moderode(A,se)
%% Modified Version of MATLAB's built-in imerode
%
% B = opterode(A,se)
%
% Details   This is a modified form of MATLAB's built-in imerode that has 
%           been optimised for 2D images with flat structural elements and 
%           was implemented for the paper 'A New Method for Ellipse 
%           Detection' by Carl J. Nelson, Philip T. G. Jackson and 
%           Boguslaw Obara in 2015.
% Inputs    A - 2D image (grayscale or BW)
%           N - structural element; must be flat
% Outputs   B - eroded 2D image (grayscale or BW)
%
% Examples:
% A = rand(20);%Starting image - random
% se = strel('disk',4,0);%Structural element - disk (N.B. flat)
% B = opterode(A,se), returns grayscale image B same size as A
%
% License   See included <a href="./LICENSE/">file</a> or visit
%           <a href="https://github.com/ChasNelson1990/...
%           	Ellipse-Detection-by-Hilbert-Edge-Detection-and-Ranging/">
%			The GitHub Repository</a>

%% Pad with Zeros (c.f. MATLAB's default of Inf)
[m,n] = size(A);
[sm,sn] = size(se);
B = zeros(m+(2*sm),n+2*sn);
B(sm+(1:m),sn+(1:n)) = A;

%% Apply Erosion
height = zeros([sm,sn]);
if islogical(A)
    B = logical(B);
    B = images.internal.morphmex('erode_binary_twod', B, se, height, -1);
else
    B = images.internal.morphmex('erode_gray_flat', B, se, height, -1);
end

%% Remove Padding
B = B(sm+(1:m),sn+(1:n));

end
