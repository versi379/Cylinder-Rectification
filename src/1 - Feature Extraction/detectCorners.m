addpath(genpath('Functions'));

I = imread('../../img/PalazzoTe.jpg');
I = imrotate(I,270);
I = rgb2gray(I);
I = im2double(I);

% find corners -- Harris Algorithm
[corner_x, corner_y] = findCorners(I, 3, 20);

% plot corners on image
figure(1);
imshow(I);
title('Corner Detection');
hold on
plot(corner_y, corner_x, append('yellow', '+')); % switched xy due to coord. system
hold off

% save image with detected corners
saveas(1, '../../img/1 - Feature Extraction/PalazzoTe_corners.png');
