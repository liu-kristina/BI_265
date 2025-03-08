% DCE Analyses
%% 1 - Setup Constants, Filenames
close all;
clear all;
%maxenhance = xxx;  
%basetime = xx;
Slice = 10;  % Select a slice to work with, e.g. 10

%Defining directories
DCE_dir = 'data/E317/1009/';
rootname = 'E317S1009I';

%Collect the name of all the images
%cd DCE_dir ;
cd 'data/E317/1009';
i=1;
ImName = strcat(rootname,num2str(i),'.DCM');
tempres = 10.471;   %temporal resolution in sec (based on data, but use as constant here)

%% - 2 Get Scan Parameters

%Extract some variables from info of the first image
DCEinfo = dicominfo(ImName);   
numfiles = DCEinfo.ImagesInAcquisition;
numtimepts = DCEinfo.NumberOfTemporalPositions;
numslices = numfiles / numtimepts;

%% 3-Load the 4D data - time, slice, image(x,y) - Note: this will take a long time

% Open files - data is stored file1=slice1, time 1; file2 = slice2, time 1 ...
for t = 0:(numtimepts-1)
    for i = 1:numslices 
        filenumber = t*numslices + i;
        ImName = strcat(rootname,num2str(filenumber),'.DCM');
        DCMinfo = dicominfo(ImName);
        ImDCE(t+1,i,:,:) = dicomread(DCMinfo);
    end
end

%% 4 - Calculate baseline intensity, Make Enhancement Images (normalize to baseline), max slope and max peak maps
%Can just use data from 1 slice
RawSl = squeeze(ImDCE(:,Slice,:,:));

% Initialize accumulation variables and max variables to 0

baseline = zeros(1,256,256);
enhancement = zeros(30,256,256);
slope = zeros(30,256,256);

% Loop over all slices  --- or analyze just for your selected slice

% Find average of baseline points
baseline = mean(RawSl(1:5,:,:));

% Loop over all timepoints
% Calculate Enhancement at each timepoint

enhancement = 100.*(double(RawSl(:,:,:))./baseline);

% Calculate the Slope at each timepoint
for i = 2:30
    slope(i,:,:) = enhancement(i,:,:)-enhancement(i-1,:,:);
end
    
% Calculate the max peak map for the slice
max_peak_map = max(enhancement);
max_peak_map = squeeze(max_peak_map(:,:,:));
% Calculate the max slope map for the slice
max_slope_map = max(slope);
max_slope_map = squeeze(max_slope_map(:,:,:));
% See help for max -- max(3Darray) --> finds the max value across the first dimention and saves a 2D array


%% 5 - Display Baseline average image, Enhancement Image at timepoint 30, and Parameter Maps for the selected slice
close all;
baseline = squeeze(baseline(:,:,:));
figure()
imagesc(baseline); colormap gray
title('Baseline Average')
enhancement30 = squeeze(enhancement(30,:,:));
figure()
imagesc(enhancement30); colormap gray
title('Enhancement at timepoint 30')
figure()
imagesc(max_slope_map, [0 100]); colormap gray
title('Max Slope Map for Slice 10')
figure()
imagesc(max_peak_map); colormap gray
title('Max Peak Map for Slice 10')

%% 6 - Create Right Peripheral Zone ROI & Left Peripheral Zone ROI -- Compare suspicious region to contralateral side
close all;
%Select slice of interest to draw Right and Left Peripheral Zone ROIs

%Manually Draw Patient's Right Peripheral Zone (Image Left)

%Display image
figure('name','Max Slope')
imagesc(max_slope_map, [0 50]);
title('Manually Draw ROI');
colormap gray
%Code to draw an ROI
freehandroi=imfreehand(gca);
enhancement_mask_right=createMask(freehandroi);
%%
%Manually Draw Patient's Left Peripheral Zone (Image Right)
figure('name','Max Slope')
imagesc(max_slope_map,[0 50]);
title('Manually Draw ROI');
colormap gray
%start timing
tic
%Code to draw an ROI
freehandroi=imfreehand(gca);
enhancement_mask_left=createMask(freehandroi);

%-------------------
%% 7-Calculate & Display ROI stats
% Calculate mean, std, median for Left ROI and Right ROI on Parameter Maps
max_peak_map(isnan(max_peak_map))=0;
max_enhancement_map(isnan(max_peak_map))=0;

c = 0;

for i = 1:256
    for j = 1:256
        if enhancement_mask_right(i,j) == 1
            c = c+1;
            right_roi_peak(c)=double(enhancement_mask_right(i,j))...
                .*max_peak_map(i,j);
        end
    end
end

c = 0;
for i = 1:256
    for j = 1:256
        if enhancement_mask_left(i,j) == 1
            c = c+1;
            left_roi_peak(c)=double(enhancement_mask_left(i,j))...
                .*max_peak_map(i,j);
        end
    end
end
c = 0;

for i = 1:256
    for j = 1:256
        if enhancement_mask_right(i,j) == 1
            c = c+1;
            right_roi_slope(c)=double(enhancement_mask_right(i,j))...
                .*max_slope_map(i,j);
        end
    end
end
c = 0;

for i = 1:256
    for j = 1:256
        if enhancement_mask_left(i,j) == 1
            c = c+1;
            left_roi_slope(c)=double(enhancement_mask_left(i,j))...
                .*max_slope_map(i,j);
        end
    end
end
%%

mean_right_peak=nanmean(mean(right_roi_peak))
mean_left_peak= nanmean(mean(left_roi_peak))
median_right_peak = nanmedian(median(right_roi_peak))
median_left_peak = nanmedian(median(left_roi_peak))
std_right_peak = nanstd(std(right_roi_peak))
std_left_peak =nanstd(std(left_roi_peak))

mean_right_slope=nanmean(mean(right_roi_slope))
mean_left_slope= nanmean(mean(left_roi_slope))
median_right_slope = nanmedian(median(right_roi_slope))
median_left_slope = nanmedian(median(left_roi_slope))
std_right_slope = nanstd(std(right_roi_slope))
std_left_slope =nanstd(std(left_roi_slope))

%% 8- Plot uptake curves vs. time in seconds

%Set x-axis values in sec
t = tempres.*[1:30];
%Calculate mean ROI vs. time
mean_roi_right= zeros(1,30);

%%
mean_roi_right=zeros(1,30);
mean_roi_left = zeros(1,30);
for i = 1:30
    enhancement_2d=squeeze(enhancement(i,:,:));
    overlay_right= enhancement_2d(enhancement_mask_right);
    overlay_left= enhancement_2d(enhancement_mask_left);
    mean_roi_right(:,i)=mean(overlay_right);
    mean_roi_left(:,i)=mean(overlay_left);
end

%%
close all
% mean enhancement on y axis, time on x axis

%Plot both curves - see plot, hold on, legend 
 
plot(t, mean_roi_right,'or')
hold on
plot(t, mean_roi_left,'ob')
legend('Right ROI','Left ROI')

%% 9 - Calculate Washout Slope for the 2 ROIs as a linear fit to the last half of the timepoints

%create x vector
xvec  = tempres.*[16:30];
%create y vector of the same size
yvec_right = mean_roi_right(16:30);
yvec_left = mean_roi_left(16:30);
p_right = polyfit(xvec, yvec_right,1);  %first order polynomial: y = mx + b format
rwo_slope = p_right(1);
p_left = polyfit(xvec, yvec_left,1);  %first order polynomial: y = mx + b format
lwo_slope = p_left(1);

% Put in units of %baseline / minute

rwo_slope = rwo_slope*60
lwo_slope = lwo_slope*60


%---------------------------------------------------------------
%% - Reset directory out of data folder
cd ../../..
%% - Answer Questions from Lab Handout
% 1. Cancer is most likely on patient's left. The left ROI had a steep
% upslope and a washout while the right enhancement increased more slowly 
% and kept increasing.
% 2. High grade, because the steeper slope means more permeability and
% higher peak enhancement means more washout of extracellular extracellular 
% space.
% 3. a) Not easy from baseline, image doe not have much contrast
% difference. b) Not from enhancement at timepoint 30 either. c) Not max
% peak map. d) Max slope map has contrast difference. e) Very clear from
% the curves.
% 4. Getting the right brightness and contrast is challenging. You also
% need to know the appropriate T1 and T2 to get that contrast in DCE in an
% appropriate time frame.
% 5. Getting the right contrast is crucial for viewing the cancer. 
% Weighting the max slope image with the T2 to get even better contrast.
% Measurements could be more accurate with better selection of region pf
% interest using that contrast.


