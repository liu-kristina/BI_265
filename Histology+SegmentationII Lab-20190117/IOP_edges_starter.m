%%Edge Detect and Region Fill to Segment SAT
%
clear all; close all; clc
%constants
%uncomment and set an absolute threshold for edge detection if use this
%method:
%edgedetectthreshold = 8;   %absolute threshold for edge detection
current_dir = pwd;
current_file = '../MRFat/data/suc047_4_S10';
VATslice = 30;  %15 if counting from first slice, 30 if counting from superior end

%Load images
%input = readidf_file_pc(current_file, current_dir);
%See help for options, esp for different operating systems
%image data will be stored in input.img

Inr = read_idf_image_pc(strcat(current_file, '_In'), current_dir, 0);
Outr = read_idf_image_pc(strcat(current_file, '_Out'), current_dir, 0);

%record image size parameters
img_size = size(Inr.img(:,:,:));
num_slices = img_size(3);

%Rotate and Flip the images to be in matlab-style format
for j=1:num_slices
    In(:,:,j) = flipud(imrotate(Inr.img(:,:,j),90));
    Out(:,:,j) = flipud(imrotate(Outr.img(:,:,j),90));
end


%% Calculate Fat Fraction Map, based on In-phase and Out-of-phase images, in [%]
%Calculate FF map, in [%]
%>>>>>>>>>>>>>>


%<<<<<<<<<<<<

%Remove NaN, Inf, -Inf from the map (here called FF, change if needed)
%see help isnan and help isinf
for i = 1:img_size(1)
   for j=1:img_size(2)
       for k = 1:img_size(3)

%>>>>>>>>>>>

%<<<<<<<<<<<
       
       end
   end
end

%% Edge Detect FF slice of interest

%Use Slice of Interest
%>>>>>>>>>>

%<<<<<<<<<<


%Edge detect, can use EdgeMap = edge(Image);
%>>>>>>>>>>>>

%see help edge
%The default uses 'sobel' and find a threshold. This threshold will likely
%be poor.  You can modify the automatic threshold or specify one.

%One method to find a threshold based on the image data
%[~, threshold] = edge(FF30, 'sobel')
%fudgeFactor = xxx;    %make a threshold relative to the automatically
%found one
%EdgeMap = edge(FF30,'sobel', threshold * fudgeFactor);

%Or use a specific threshold (constant listed at the top of the m-file.
%
%Edge detect FFmap (with no NaN & no Inf)
%EdgeMap = edge(FF30,'sobel', edgedetectthreshold);

%<<<<<<<<<<<<


%Display FF map and Edge detected mask
%>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<

%Define structure elements -- basic shapes to identify in the image
%e.g. are lines, - can make vertical, horizontal, 2 diagonals (+45deg and -45deg)
%set the size of each element (can vary) and orientation.
%Can try other shapes, i.e. 'disk' (filled circle)

%>>>>>>>>>>>>>>

help strel
%format:  se90 = strel('line', size, degrees); 
%e.g.  se90 = strel('line', 4, 90); 
%
%
%This element is a disk = filled circle:
%sedisk = strel('disk',size);

%<<<<<<<<<<<<

%Dilate the FF slice and the edgemap using struct elements:
%>>>>>>>>>>>>>>>>>
help imdilate
%EdgeMapdil = imdilate(EdgeMap, [list of struct elements ]);

%<<<<<<<<<<<<<<<<<<<<

%Display the dilated images
%>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<

%%Erode the images and display
help imerode;
%>>>>>>>>>>>>>>>



%<<<<<<<<<<<<<<<<<<

%% Region Grow & Display the 2 eroded, dilated images
%format: [polygonofROI, MaskofROI] =
%regionGrowing(Im,[InitialPositionx,y],1,threshold);  %leave the FF type as
%'1' (3rd input)
%borderthresh = xx;  %set a value for the cutoff. This will be different
%for the 2 maps
%change "Image" to your mask name
[SATpoly,SATmask] = regionGrowing(Image,[175,75],1,borderthresh);



%%


