function [data] = hedar (im,maxLength,th,pruning)
%% hedar, Ellipse Detection by Hilbert-Edge Detection and Ranging
%
% data = hedar (im,maxLength,th)
% data = hedar (im,maxLength,th,pruning)
%
% Details   This function contains the an ellipse detection algorithm
%			as implemented for 'Combining Mathematical Morphology
%           and the Hilbert Transform for Fully Automatic Nuclei
%           Detection in Fluorescence Microscopy' by Carl J. Nelson, 
%           Philip T. G. Jackson and Boguslaw Obara (submitted LNCS).
%			This is an advanced method that utilised the Hilbert
%			transform to find steps in the erosion signals.
%
% Inputs    im - a 2D binary or greyscale image
%           maxLength - maximum major axis length to search for (in
%           pixels; default is the length of the smallest image dimension;
%           th - the threshold (within 0 and 1) over which any step in
%           signal must be to contribute to the ellipse detection stage
%
% Outputs   data - a matrix of ellipses. Each row contains five elements:
%           the center of the ellipse, its major and minor axes and
%           orientation of the major axis.
%
% Examples:
% data = hedar (im,20,0.2), runs hedar on im looking for ellipses of maximum
% size 20 pixels and considers all signals with a step of >=20% the maximum
% to be contributing.
%
% Copyright 2019 Carl J. Nelson, Durham University, UK
%
% License   See included <a href="./LICENSE/">file</a> or visit
%           <a href="https://github.com/ChasNelson1990/...
%           	Ellipse-Detection-by-Hilbert-Edge-Detection-and-Ranging/">
%			The GitHub Repository</a>

%% Inputs
if nargin<4; pruning=true; end
if nargin<3; th = 0.2; end
if nargin<2 || isempty(maxLength) || max(maxLength==0); maxLength = min(size(im)); end

%% Set Up
im = imadjust(im);  % should this be needed?!
[m,n] = size(im);
lengthSteps = 1:maxLength;
angularSteps = linspace(0,180,min(2*maxLength,180));  % use appropriate number of angular steps but never exceed resolution of 1 degree
angularSteps(end) = [];  % 0 degrees == 180 degrees due to symettry
timing = ones(4,1);

%% Granulometric Signals
granulometricSignals = zeros(m, n, length(lengthSteps),length(angularSteps));
for lStep = 1:length(lengthSteps)  % for each length
    for aStep = 1:length(angularSteps)  % for each angle
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
        granulometricSignals(:,:,lStep,aStep) = imerode(im,se);
    end
end
clear lStep aStep lStepNumber aStepNumber se t

%% Find Drop in Signals using Hilbert
signals = reshape(granulometricSignals,[],size(granulometricSignals,3),size(granulometricSignals,4));
hilbdrop=zeros(m*n,size(granulometricSignals,4));
gap=zeros(m*n,size(granulometricSignals,4));
for ti = 1:size(signals,3)
    signal = squeeze(signals(:,:,ti))';%all signals in this orientation
    h = hilbert([signal(1,:);signal;signal(end,:)]);%Analytic Signal, note the boundary padding
    hi = imag(h(2:end-1,:));%Hilbert Transform (imaginary component of AS), note removes boundary padding
    [~,hmi] = max(hi);%position of maximum
    hmi(hi(1,:)>=0)=1;%catch for pixels caught between two ellipses
    mask = repmat((1:size(granulometricSignals,3))',1,size(signal,2));
    mask2 = repmat(hmi,size(signal,1),1);
    mask = (mask<mask2);%mask of signals below their respective maximum index
    high = sum(signal.*mask)./hmi;%average value of first part of step
    low = sum(signal.*(~mask))./-(hmi-size(signals,3));%average value of lower part of step
    low(hmi==size(signals,3))=0;%catch for infite signals
    gap(:,ti) = high-low;%size of step gap
    hilbdrop(:,ti) = squeeze(hmi)';%step position
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
minorOrientation = mod(majorOrientation-floor(length(angularSteps)/2),length(angularSteps))+1;
[y,x] = ndgrid(1:m,1:n);
idr = sub2ind(size(hilbdrop),y(:),x(:),minorOrientation(:));
minorAxis = hilbdrop(idr);
minorAxis = reshape(minorAxis,size(majorAxis));
clear minorOrientation hilbdrop x y idr idx X t granulometricSignals

%% Background Check
backgroundM = (majorAxis==1); backgroundm = (minorAxis==1);
backgroundG = gap<=(max(gap(:))*th);%(gap<=(min(gap(:))+((max(gap(:))-min(gap(:)))*th)));
background = backgroundG | backgroundM | backgroundm;
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
%regionalMaxima = logical(regionalPeaks - imreconstruct(regionalPeaks-1,regionalPeaks));
%centroids = regionprops(regionalMaxima,'Centroid');
%centroids = reshape(cell2mat(struct2cell(centroids))',2,[])';
regionalMaxima = imregionalmax(regionalPeaks);
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
data(:,1:2) = data(:,1:2)+1;  % hack
clear idc centroids majorAxis minorAxis majorOrientation lengthSteps angularSteps

%% Prune Duplicates by Overlap
dists = pdist2(data(:,1:2),data(:,1:2));
if pruning==true
    marked = zeros(size(data,1),1);
    for i=1:size(data,1)
        ellipsei = ellipse2([m,n],data(i,1:2),data(i,3),data(i,4),data(i,5));%2.4
        ellipsei = logical(ellipsei);
        if sum(ellipsei)<pi
            marked(i) = 1;
            continue;
        end
        for j=i+1:size(data,1)
            if marked(i)==0 && marked(j)==0 && dists(i,j)<=2*maxLength
                ellipsej = logical(ellipse2([m,n],data(j,1:2),data(j,3),data(j,4),data(j,5)));
                if sum(ellipsej)<pi
                    marked(j) = 1;
                    continue;
                end
                common = sum(ellipsei(:) & ellipsej(:));
                score1 = common/sum(ellipsei(:));
                score2 = common/sum(ellipsej(:));
                if max(score1,score2)>0.5
                    if score1>score2
                        marked(j) = 1;
                    else
                        marked(i) = 1;
                    end
                end
            end
        end
    end
    data(logical(marked),:) = [];
    clear i j ellipsei ellipsej marked common score1 score2
end

end
