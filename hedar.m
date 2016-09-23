function data = hedar (im,maxLength,th,resLength,resAngular)
%% hedar, Ellipse Detection by Hilbert-Edge Detection and Ranging
%
% data = hedar (im,maxLength,th)
% data = hedar (im,maxLength,th,resLength,resAngular)
% data = hedar (im,[minLength,maxLength],...)
%
% Details   This function contains the an ellipse detection algorithm
%			as implemented for 'Ellipse Detection by Hilbert-Edge
%			Detection and Ranging' by Carl J. Nelson, Philip T. G.
%			Jackson and Boguslaw Obara (submitted for publication).
%			This is an advanced method that utilised the Hilbert
%			transform to find steps in the erosion signals.
%
% Inputs    im - a 2D binary or greyscale image
%           maxLength - maximum major axis length to search for (in
%           pixels; default is the length of the smallest image dimension;
%           if this argument is a vector then minLength = maxLength(1),
%           maxLength=maxLength(2), thus a minimum axis length can be set)
%           th - the threshold (within 0 and 1) over which any step in
%           signal must be to contribute to the ellipse detection stage
%           resLength - the resolution to search for axis lengths between
%           minLength and maxLength (default is 1 pixel)
%           resAngular - the resolution to search for axis orientation
%           between 0 and 180 degrees (default is 1 degree)
%
% Outputs   data - a matrix of ellipses. Each row contains five elements:
%           the center of the ellipse, its major and minor axes and
%           orientation of the major axis.
%
% Examples:
% data = hedar (im,20,5), runs hedar on im looking for ellipses of maximum
% size 20 pixels and resolution of +/- 5 pixels and with angular
% resolution of 45 degrees (the default), considers all signals with a step
% of >=20% (default) the maximum to be contributing
% data = hedar (im,[10,20],[],5,0), runs hedar on im looking for ellipses of
% axis size between 10 and 20 pixels and angular resolution of 5 degrees,
% considers all signals with a step of >=0% the maximum to be contributing,
% i.e. all steps
%
% Copyright 2016 Carl J. Nelson, Durham University, UK
%
% License   See included <a href="./LICENSE/">file</a> or visit
%           <a href="https://github.com/ChasNelson1990/...
%           	Ellipse-Detection-by-Hilbert-Edge-Detection-and-Ranging/">
%			The GitHub Repository</a>

%% Inputs
if nargin<5 || isempty(resAngular); resAngular = 1; end
if nargin<4 || isempty(resLength); resLength = 1; end
if nargin<3; th = 0.2; end
if nargin<2 || isempty(maxLength) || max(maxLength==0); maxLength = min(size(im)); end
if length(maxLength)>1
    minLength = maxLength(1);
    maxLength = maxLength(2);
else
    minLength = 1;
end

%% Set Up
[m,n] = size(im);
lStepNumber = ceil((maxLength-minLength+1)/resLength);
lengthSteps = linspace(minLength,maxLength,lStepNumber);
aStepNumber = ceil(180/resAngular);
angularSteps = linspace(0,180,aStepNumber+1);
angularSteps(end) = [];

%% Granulometric Signals
granulometricSignals = zeros(m, n, lStepNumber,aStepNumber);
for lStep = 1:length(lengthSteps)
    for aStep = 1:length(angularSteps)
        %Structuring Element
        r = floor(lengthSteps(lStep)/2); a = angularSteps(aStep);
        p1 = [-r*sind(a),r*cosd(a)];
        p2 = [r*sind(a),-r*cosd(a)];
        pmin = min([p1;p2]);
        p1 = round(p1-pmin)+1; p2 = round(p2-pmin)+1;
        pmax = max([p1;p2]);
        se = false(pmax); se(p1(1),p1(2)) = true; se(p2(1),p2(2)) = true;
        clear r a p1 p2 pmin pmax
        %Erosion
        granulometricSignals(:,:,lStep,aStep) = moderode(im,se);
    end
end
clear lStep aStep lStepNumber aStepNumber se

%% Find Drop in Signals using Hilbert
signals = reshape(granulometricSignals,[],size(granulometricSignals,3),size(granulometricSignals,4));
hilbdrop=zeros(m*n,size(granulometricSignals,4));
gap=zeros(m*n,size(granulometricSignals,4));
for ti = 1:size(signals,3)
    signal = squeeze(signals(:,:,ti))';
    h = hilbert([signal(1,:);signal;signal(end,:)]);
    h = h(2:end-1,:);
    hi = imag(h);
    [~,hmi] = max(hi);
    hmi(hi(1,:)>=0)=1;%catch for pixels caught between two ellipses
    mask = repmat((1:size(granulometricSignals,3))',1,size(signal,2));
    mask2 = repmat(hmi,size(signal,1),1);
    mask = (mask<mask2);
    high = sum(signal.*mask)./hmi;
    low = sum(signal.*(~mask))./-(hmi-size(signals,3));
    gap(:,ti) = high-low;
    hilbdrop(:,ti) = squeeze(hmi)';
    clear signal ti h hi hmi mask mask2 high low
end
hilbdrop = reshape(hilbdrop,size(granulometricSignals,1),size(granulometricSignals,2),size(granulometricSignals,4));
clear signals

%% Major & Minor Axes
[majorAxis,majorOrientation] = max(hilbdrop,[],3);
X = 1:size(gap,1);
idx = sub2ind(size(gap),X',majorOrientation(:));
gap = gap(idx);
gap = reshape(gap,size(granulometricSignals,1),size(granulometricSignals,2));
minorOrientation = mod(majorOrientation-((90+resAngular)/resAngular),180/resAngular)+1;
[y,x] = ndgrid(1:m,1:n);
idr = sub2ind(size(hilbdrop),y(:),x(:),minorOrientation(:));
minorAxis = hilbdrop(idr);
minorAxis = reshape(minorAxis,size(majorAxis));
clear minorOrientation hilbdrop x y idr idx X granulometricSignals

%% Background Check
backgroundM = (majorAxis==1); backgroundm = (minorAxis==1);
backgroundG = gap<=(max(gap(:))*th);%(gap<=(min(gap(:))+((max(gap(:))-min(gap(:)))*th)));
background = backgroundM | backgroundm | backgroundG;
majorAxis = majorAxis-1; minorAxis = minorAxis-1;
majorOrientation = majorOrientation-1;
majorAxis(background) = 1; majorAxis(majorAxis==0) = 1;
majorOrientation(background) = 1; majorOrientation(majorOrientation==0) = 1;
minorAxis(background) = 1; minorAxis(minorAxis==0) = 1;
clear backgroundM backgroundm backgroundG background

%% Centroids
minorAxisMasked =  minorAxis.* (minorAxis>1);
majorAxisMasked =  majorAxis.* (majorAxis>1);
regionalPeaks = majorAxisMasked.*minorAxisMasked;
regionalPeaks = medfilt2(regionalPeaks, [3,3]);
regionalMaxima = logical(regionalPeaks - imreconstruct(regionalPeaks-1,regionalPeaks));
centroids = regionprops(regionalMaxima,'Centroid');
centroids = reshape(cell2mat(struct2cell(centroids))',2,[])';
clear minorAxisMasked majorAxisMasked regionalPeaks regionalMaxima

%% Compile Data
data = zeros(size(centroids,1),5);
data(:,1:2) = round(centroids);
idc = sub2ind([m,n],data(:,2),data(:,1));
data(:,3) = lengthSteps(majorAxis(idc));
data(:,4) = lengthSteps(minorAxis(idc));
data(:,5) = angularSteps(majorOrientation(idc));
clear idc centroids majorAxis minorAxis majorOrientation lengthSteps angularSteps

%% Prune Duplicates by Overlap
marked = zeros(size(data,1),1);
for i=1:size(data,1)
    ellipsei = logical(ellipse2([m,n],data(i,1:2),data(i,3),data(i,4),data(i,5)));
    if sum(ellipsei)<pi*minLength^2
        marked(i) = 1;
        continue;
    end
    for j=i+1:size(data,1)
        if marked(i)==0 && marked(j)==0 && pdist2(data(i,1:2),data(j,1:2))<=2*maxLength
            ellipsej = logical(ellipse2([m,n],data(j,1:2),data(j,3),data(j,4),data(j,5)));
            common = sum(ellipsei(:) & ellipsej(:));
            score1 = common/sum(ellipsei(:));
            score2 = common/sum(ellipsej(:));
            if max(score1,score2)>0.5
                if score1>score2
                    marked(i) = 1;
                else
                    marked(j) = 1;
                end
            end
        end
    end
end
data(logical(marked),:) = [];
clear i j ellipsei ellipsej marked common score1 score2

end
