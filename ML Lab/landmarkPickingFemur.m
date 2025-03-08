function landmarkPickingFemur (inputFileName, p)

% Script to pick landmarks on 2D images
% by Serena Bonaretti, UCSF

% Input:
% - inputFileName: file name of the image where to pick the landmarks
% - p: number of landmarks to pick
% Output:
% - *coords.txt: vector that contains the coordinates of the picked landmarks. 
%                The first half of the vector contains the x, the second half contains the y. 

% close all opened figure
close all

% load image
img = imread (inputFileName);
figure, imagesc (img), axis image

% pick the landmarks
for i =1:p
    [x(i),y(i)] = ginput(1);
    h1 = text(x(i),y(i),'*', 'HorizontalAlignment','center', 'Color', [0 0 1], 'FontSize',16);
    h2 = text(x(i),y(i),num2str(i), 'HorizontalAlignment','center', 'Color', [1 0 0], 'FontSize',14);
end

% write coordinates to txt file
[pathstr,name,ext] = fileparts(inputFileName); 
outputFileName = [pathstr '/' name '_coords.txt'];
[fid msg] = fopen(outputFileName, 'w');

for i = 1:length(x)
    fprintf(fid, '%.1f\n', x(i));
end
for i = 1:length(y)
    fprintf(fid, '%.1f\n', y(i));
end
      
fclose(fid);

% double check saved coordinates
figure, imagesc (img), axis image
coords = load(outputFileName);
hold on; plot(coords(1:length(coords)/2), coords(length(coords)/2+1:end),'*')