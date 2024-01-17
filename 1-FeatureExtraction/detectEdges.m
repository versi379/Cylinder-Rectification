I = imread('../img/PalazzoTe.jpg');
I = imrotate(I,270);
I = rgb2gray(I);
I = im2double(I);

% Otsu Method -- select threshold
threshold = graythresh(I)

% Canny Algorithm -- edge detection
edges = edge(I, 'canny', [0.025 0.05]);

% show image with detected edges
figure(1),
imshow(edges);
title('Edge Detection');

% save image with detected edges
saveas(edges, '../img/Extracted Features/PalazzoTe_edges.png');
