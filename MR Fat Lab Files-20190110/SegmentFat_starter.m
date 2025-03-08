clear all; close all; clc
%Part I - V
%Load images
current_dir = pwd;

%input = readidf_file_pc(current_file, current_dir);
%See help for options, esp for different operating systems
%image data will be stored in <input>.img

%Pre Diet 
current_file = 'data/suc047_4_S10';
VATslice = 30;  %15 if counting from first slice, 30 if counting from superior end

Inr = read_idf_image_pc(strcat(current_file, '_In'), current_dir, 0);
Outr = read_idf_image_pc(strcat(current_file, '_Out'), current_dir, 0);
Fatr = read_idf_image_pc(strcat(current_file, '_F'), current_dir, 0);
Wr = read_idf_image_pc(strcat(current_file, '_W'), current_dir, 0);

%Post Diet
%current_file = 'data/suc047_2_S21';
%VATslice = 27;   %18 if counting from first slice (inferior), 27 if counting
%from superior end

%record image size parameters 
img_size = size(Inr.img(:,:,:))
%FYI
%img_rl = img_size(1);
%img_ap = img_size(2);
num_slices = img_size(3);
%other methods of collecting information: im.data, im.img

%Rotate and Flip the images to be in matlab-style format
for j=1:num_slices
    F(:,:,j) = flipud(imrotate(Fatr.img(:,:,j),90));
    W(:,:,j) = flipud(imrotate(Wr.img(:,:,j),90));
    In(:,:,j) = flipud(imrotate(Inr.img(:,:,j),90));
    Out(:,:,j) = flipud(imrotate(Outr.img(:,:,j),90));
end


%Display the image by slice - check read correctly
figure('Name', strcat('Original Image: ', current_file))
for j=1:img_size(3)
    imagesc(F(:,:,j))
    axis image
    title({'Fat: Slice',j});
    colormap gray
    pause(0.1)
end

%%
%Part I - 2: Calculate Fat Fraction Map in [%]

%Create Body/Good SNR Mask


%Threshold an image set
%example:
Inmask = In>100;

%Apply a 2D filter to the mask (see help medfilt2)


%Evaluate filter -> Display, confirm threshold and filter choice

%Apply this Body/SNR Mask to the FF map


%Select Slice of Interest
%example:
FFi=FF(:,:,VATslice);


%%
%Part II-6 - Display FF, Body/SNR mask, mFF

%Manually Draw Visceral Fat ROI
%Display image
figure('name','FF')
imagesc(FFi);
title('Manually Draw ROI of VAT');
colormap gray
%start timing
tic
%Code to draw an ROI
freehandroi=imfreehand(gca);
manVATmask=createMask(freehandroi);

%end timing
manualVATtime = toc

%Part II -2 - Save image w/ drawn ROI

%Characterize VAT ROI - histogram, size in pixels, mean, median, std
help histogram

%eg:
manVATmean = mean(FFi(manVATmask)),



%%
%Part III  - Semiautomatically Generate VAT ROI - 

%1-Manually Draw Visceral ROI

%Create a mask

%Create a histogram of the FF in the viscera

%%
%2-Automatically generate VAT mask and VAT-FF for:
%   1-Use Fat Fraction Map, 2 - Use masked Fat Fraction Map

%Threshold the FF map to > 50

%Generate a VAT ROI Fat Fraction (masked by viscera and by FF > 50)

%Create a histogram of the VAT ROI Fat Fraction values

%Make a Figure of VAT ROI FF image

%--Repeat for the masked Fat Fraction Map


%%
% Part IV - Determine Accuracy, Differences
%A - Compare VAT ROI on masked FF vs. manual VAT ROI
%B - Compare VAT ROI on masked FF to the VAT ROI on FF (not masked)

%Calculate Dice Coefficient & the %Difference in Size

%Display FF image, 2 masks, and overlap of the 2 masks

%Calculate & Report stats on the VAT ROIs
%Summary of Results
disp(sprintf('Parameter        Manual VAT        Semiauto VAT'))


%%
%Part V - Compare Post VAT vs. Pre VAT

