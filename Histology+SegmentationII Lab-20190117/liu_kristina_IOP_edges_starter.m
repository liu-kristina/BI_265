%%Edge Detect and Region Fill to Segment SAT
%
clear all; close all; clc
%constants
%uncomment and set an absolute threshold for edge detection if use this
%method:
%edgedetectthreshold = 8;   %absolute threshold for edge detection
current_dir = pwd;
current_file = 'data/suc047_4_S10';
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
close all;
FF = (In(:,:,VATslice)-Out(:,:,VATslice))./(2*In(:,:,VATslice));
FF(isnan(FF))=0;
FF(isinf(FF))=0;
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
close all;
%Use Slice of Interest
%>>>>>>>>>>
close all;
MaskIn=zeros(img_size(1));
for i = 1:img_size(1);
    for j = 1:img_size(2);
        if In(i,j,30)>100
            MaskIn(i,j)=1;
        else 
            MaskIn(i,j)=0;
        end
    end
end
MaskIn=medfilt2(MaskIn);
filtFF = MaskIn.*FF;
figure('Name', 'Filtered FF')
imagesc(filtFF); colormap gray


%%
%Edge detect, can use EdgeMap = edge(Image);
%>>>>>>>>>>>>
[~, threshold] = edge(filtFF, 'sobel');
EdgeMap = edge(filtFF, 'sobel', 0.08);

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

close all;
%Display FF map and Edge detected mask
%>>>>>>>>>>>>>>>
figure('Name', 'FF map and Edge mask')
subplot(2,1,1)
title('Fat Fraction')
imagesc(filtFF); colormap gray
subplot(2,1,2)
imagesc(EdgeMap); colormap gray
title('Edge map')


%<<<<<<<<<<<<<<<<
%%
close all;
%Define structure elements -- basic shapes to identify in the image
%e.g. are lines, - can make vertical, horizontal, 2 diagonals (+45deg and -45deg)
%set the size of each element (can vary) and orientation.
%Can try other shapes, i.e. 'disk' (filled circle)

%>>>>>>>>>>>>>>

% help strel
%format:  se90 = strel('line', size, degrees); 
%e.g.  se90 = strel('line', 4, 90); 
%
%
%This element is a disk = filled circle:
%sedisk = strel('disk',size);

%<<<<<<<<<<<<

%Dilate the FF slice and the edgemap using struct elements:
%>>>>>>>>>>>>>>>>>
% help imdilate
se_dil = strel('disk', 2);
EdgeMapdil = imdilate(EdgeMap, se_dil);
se_dil_FF = strel('disk', 2);
FFdil = imdilate(filtFF, se_dil_FF);
%<<<<<<<<<<<<<<<<<<<<

%Display the dilated images
%>>>>>>>>>>>>>>>>>>
figure('Name', 'Dilated Edge Map')
imagesc(EdgeMapdil); colormap gray
figure('Name', 'Dilated FF Map')
imagesc(FFdil); colormap gray
%<<<<<<<<<<<<<

%%Erode the images and display
% help imerode;
%>>>>>>>>>>>>>>>

se_er = strel('disk', 1);
EdgeMaper = imerode(EdgeMapdil, se_er);
se_er_FF = strel('disk', 1);
FFer = imerode(FFdil, se_er_FF);

figure('Name', 'Eroded Edge Map')
imagesc(EdgeMaper); colormap gray
title('Edge Map After Erosion')
figure('Name', 'Eroded Fat Fraction Map')
imagesc(FFer); colormap gray
title('Fat Fraction After Erosion')


%<<<<<<<<<<<<<<<<<<

%% Region Grow & Display the 2 eroded, dilated images
close all;
%format: [polygonofROI, MaskofROI] =
%regionGrowing(Im,[InitialPositionx,y],1,threshold);  %leave the FF type as
%'1' (3rd input)
%borderthresh = xx;  %set a value for the cutoff. This will be different
%for the 2 maps
%change "Image" to your mask name
[SATpoly_edge,SATmask_edge] = regionGrowing(EdgeMaper,[175,75],1, 0);
SATmask_edge = double(SATmask_edge);
SATmask_edge_rg= SATmask_edge.*filtFF;
figure('Name', 'SATmask Edgemap')
imagesc(SATmask_edge_rg, [0 0.2]); colormap gray
title('Edgemap Region Growing SAT')
[SATpoly_ff,SATmask_ff] = regionGrowing(FFer,[175,75],1, 0.172);
SATmask_ff = double(SATmask_ff);
SATmask_ff_rg = SATmask_ff.*filtFF;
figure('Name', 'SATmask FFmap')
imagesc(SATmask_ff_rg, [0 0.2]); colormap gray
title('Fat Fraction Region Growing SAT ')



%%
close all;
figure('Name', 'Visual Assessment of Edgemap')
subplot(2,1,1)
imagesc(SATmask_edge_rg, [0 0.2]); colormap gray
title('Edge Detected Masked Region Growing')
subplot(2,1,2)
imagesc(filtFF, [0 0.5]); colormap gray
title('Fat Fraction')

figure('Name', 'Visual Assessment of Fat Fraction')
subplot(2,1,1)
imagesc(SATmask_ff_rg, [0 0.2]); colormap gray
title('Fat Fraction Masked Region Growing')
subplot(2,1,2)
imagesc(filtFF, [0 0.5]); colormap gray
title('Fat Fraction')

%%
% 1. Using a larger dilation and erosion radius made the region growing with
% the FF map difficult. Changing the threshold by a small amount made
% much of the image appear, including the spine. It makes the intensities
% less differentiable.
% 2. The dilation and erosion creates a smooth connected line for the 
% region growing of the SAT. I used a disk of size 2 for the dilation and 1
% for the erosion.
% 3. The resulting SAT missed a part of the top and was overall thinner than
% in the original FF map.
% 4. By visual examination, the edge detected mask is more faithful to the
% subcutaneous fat in the fat fraction map.
% Since I used the same dilation and erosion, I got both maps to
% be similar shaped. However, the edge detected mask did not have the 
% "cuts" in the edge that were visible in the fat fraction masked.