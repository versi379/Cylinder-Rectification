addpath(genpath('Functions'));

I = imread('../../img/PalazzoTe.jpg');
I = imrotate(I,270);
I = rgb2gray(I);
I = im2double(I);

% find lines -- Hough Transform
lines = findLines(I);

% plot lines on image
figure(1);
imshow(I);
title('Line Detection');
hold on
for k = 1:length(lines) % for each line detected
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1), xy(:,2), 'LineWidth',1, 'Color','red');
    text(lines(k).point1(1), lines(k).point1(2), int2str(k), 'FontSize',15, 'Color','white');
end

% save image with detected lines
saveas(1, '../../img/1 - Feature Extraction/PalazzoTe_lines.png');
