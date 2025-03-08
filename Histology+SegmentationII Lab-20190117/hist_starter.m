%% Find Fat Droplets in Histology Slide
% Susan Noworolski
% 1/21/2018
% Use this code as a starter. Places you need to modify the code are marked
% with the following (and may have some lines started, but commented out):
%>>>>>>>>>>

%<<<<<<<<<<<<<

%% Initialize
clear all; close all;
Im=imread('histology_droplets.jpg');
%Im=imresize(Imin,4);
figure('name','Histology Slide');
imagesc(Im);
%Set min,max radii and sensitivities for imfindcircles:
%>>>>>>>>
minr = 19;  %minimum big droplet radius to identify
maxr = 60;  %maximum big droplet radius to identify
sens = 0.8; %sensitivity for identifying big droplets

%<<<<<<<<

%% Find Big Droplet Fat
figure('name','Histology Slide - Large Droplet Fat');
imagesc(Im);
[centersbig,radiibig,metricbig] = imfindcircles(Im,[minr,maxr], 'Sensitivity',sens);

% - Display big droplets on array
viscircles(centersbig, radiibig, 'Edgecolor', 'g');

% Calculate # big circles, mean radius, median radius, histogram big droplet radii
%>>>>>>>>>>>



%<<<<<<<<<<

%% Set up blank array to paint segmented circles to calculate areas

%>>>>>>>>>

%<<<<<<<<<


% Some circles are not fully within Im array (fall off edges) --> increase size of array to
% include complete circles
%
%>>>>>>>>>
%offset=???;
k=round(offset/2);

%Set up arrays to paint droplets
blank = zeros(size(Im,1)+offset, size(Im,2)+offset);

%% Paint Big Droplets on Array
%Loop through all circles, painting circles as '1', background has '0'
bigdrop = blank; % start with a blank slate
for i=1:size(centersbig,1)
    th = linspace(0,2*pi,3000); % Step through enough angles to get ~full circumference
    r=radiibig(i);  
    %Fill manually, issue w/ round off errors; may not get all pixels
    %Burn pixels as concentric circles, from center to radius
    for j=1:round(r)*2  
        x = j/2.0*cos(th) +centersbig(i,1);
        y = j/2.0*sin(th) +centersbig(i,2);
        for w = 1:length(y)
            row=round(x(w))+k;
            col=round(y(w))+k;
            bigdrop(col,row)=1;
        end
    end
end

%-------------------------------
figure('name','Blank image with filled circles');
imagesc(bigdrop);
colormap gray;

%% Find Small Droplet Fat
% Repeat similar to the Find Large Droplet Fat section; Use different
% radii, possibly sensitivity (Can copy code and modify with new variable
% names or make a function/s)

%>>>>>>>>>>>>>
%Find Small Droplet Fat

%Display circles

% Calculate # small circles, mean radius, histogram small droplet radii

%<<<<<<<<<<<<<<
%% Paint Small Droplets on Array
%Repeat similar to Paint Large Droplets on  a blank array (Can copy code and modify with new variable
% names or make a function/s)
%>>>>>>>>>>>>>>>>>>
%Loop through all circles, painting circles


%Display blank array with filled small droplet circles

%<<<<<<<<<<<<<<<<<

%% Calcuate Stats on segmented circles / droplets 
%You may have different variable names than listed for the small droplet
%fat variables. Update as needed.

%Zoom to cut off border
cropx = size(smdrop,1)-k;
cropy = size(smdrop,2)-k;
%Make mask arrays the same size as the original image 
%(circles = '1', background = '0')
cropsmdrop = smdrop(k:cropx, k:cropy);
cropbigdrop = bigdrop(k:cropx, k:cropy);

%Calculate area of small & large droplets in pixels and in % of the slide
%pixels
%>>>>>>>>>>

%Calculate totals: number, area, average radii, median radii, and total fat %area (of slide
%area):

%<<<<<<<<<<<<<


%% Summary of Large and Small Droplet Fat Segmentation
%Display both large and small droplets on the histology slide
figure('name','Histology Slide - Large (green) and Small (blue) Droplet Fat');
imagesc(Im);
%>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<

%Display a summary of results
disp(sprintf('Size  Number Total_Area Avg_Radius Median_Radius %%Slide'))
disp(sprintf('Large %5d  %9d   %7.3f     %7.3f     %5.1f %%',numbig(1), bigfat, avgbigradii, medianbigradii, bigperc))

%Repeat for small droplets/circles and for total circles
%>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<


%% Questions
%Answer the questions here or separately - whichever is easier.
