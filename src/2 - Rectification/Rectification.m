close all
clear
clc
addpath(genpath('Functions'));
%% READ IMAGE
img = imread("../../img/PalazzoTe.jpg");
img = imrotate(img,270);
figure(1), imshow(img);
hold on;
%% CROSS SECTION C1

[x, y] = ginput(); % select 5 points that lay on a conic in the image
scatter(x,y,40,'filled','o','MarkerFaceColor','b');

% conic equation -- cartesian coordinates
A = [x.^2 x.*y y.^2 x y ones(size(x))];
N = null(A);

% conic coefficients
cc = N(:, 1);
[a1, b1, c1, d1, e1, f1] = deal(cc(1),cc(2),cc(3),cc(4),cc(5),cc(6));
C1 = [a1 b1/2 d1/2; b1/2 c1 e1/2; d1/2 e1/2 f1]; % conic coefficients matrix
%% CROSS SECTION C2

[x, y] = ginput(); % select 5 points that lay on a conic in the image
scatter(x,y,40,'filled','o','MarkerFaceColor','b');

% conic equation -- cartesian coordinates
A = [x.^2 x.*y y.^2 x y ones(size(x))];
N = null(A);

% conic coefficients
cc = N(:, 1);
[a2, b2, c2, d2, e2, f2] = deal(cc(1),cc(2),cc(3),cc(4),cc(5),cc(6));
C2 = [a2 b2/2 d2/2; b2/2 c2 e2/2; d2/2 e2/2 f2]; % conic coefficients matrix
%% CROSS SECTION C3

[x, y] = ginput(); % select 5 points that lay on a conic in the image
scatter(x,y,40,'filled','o','MarkerFaceColor','b');

% conic equation -- cartesian coordinates
A = [x.^2 x.*y y.^2 x y ones(size(x))];
N = null(A);

% conic coefficients
cc = N(:, 1);
[a3, b3, c3, d3, e3, f3] = deal(cc(1),cc(2),cc(3),cc(4),cc(5),cc(6));
C3 = [a3 b3/2 d3/2; b3/2 c3 e3/2; d3/2 e3/2 f3]; % conic coefficients matrix

% save image with extracted conics
% saveas(1, '../../img/2 - Rectification/PalazzoTe_rec_conics.png');
%% CONIC EQUATIONS
syms 'x';
syms 'y';

% C1,C2,C3 equations -- cartesian coordinates
eq1 = a1*x^2 + b1*x*y + c1*y^2 + d1*x + e1*y + f1;
eq2 = a2*x^2 + b2*x*y + c2*y^2 + d2*x + e2*y + f2;
eq3 = a3*x^2 + b3*x*y + c3*y^2 + d3*x + e3*y + f3;

% intersect pair of conics
S12 = solve([eq1 == 0, eq2 == 0], [x,y]); % C1 int. C2
s1 = [double(S12.x(1));double(S12.y(1));1];
s2 = [double(S12.x(2));double(S12.y(2));1];
s3 = [double(S12.x(3));double(S12.y(3));1];
s4 = [double(S12.x(4));double(S12.y(4));1];
S13 = solve([eq1 == 0, eq3 == 0], [x,y]); % C1 int. C3
p1 = [double(S13.x(1));double(S13.y(1));1];
p2 = [double(S13.x(2));double(S13.y(2));1];
p3 = [double(S13.x(3));double(S13.y(3));1];
p4 = [double(S13.x(4));double(S13.y(4));1];
S23 = solve([eq2 == 0, eq3 == 0], [x,y]); % C2 int. C3
q1 = [double(S23.x(1));double(S23.y(1));1];
q2 = [double(S23.x(2));double(S23.y(2));1];
q3 = [double(S23.x(3));double(S23.y(3));1];
q4 = [double(S23.x(4));double(S23.y(4));1];
%% CIRCULAR POINTS IMAGE
II = mean([s1,p1,q1],2);
JJ = mean([s2,p2,q2],2);
%% DUAL CONIC TO CIRCULAR POINTS IMAGE
imDCCP = II*JJ' + JJ*II';
imDCCP = imDCCP./norm(imDCCP);
%% LINE AT INFINITY
l_inf = null(imDCCP);
%% AFFINE RECTIFICATION MATRIX
H = [eye(2), zeros(2,1); l_inf(:)'];
%% CONIC PROJECTIVE TRANSFORMATION
Q = inv(H)'*C1*inv(H);
Q = Q./Q(3,3);
%% CONIC GEOMETRIC INFORMATION EXTRACTION
par_geo = AtoG([Q(1,1),2*Q(1,2),Q(2,2),2*Q(1,3),2*Q(2,3),Q(3,3)]);
center = par_geo(1:2);
axes = par_geo(3:4);
a = axes(1);
b = axes(2);
alpha = par_geo(5);
%% ROTATION AND SCALING
U = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];
S = diag([1, a/b]);
%% COMPOSITION HOMOGRAPHY
K = U*S*U';
A = [K zeros(2,1); zeros(1,2), 1];
T = A*H;
tform = projective2d(T');
J = imwarp(img, tform);
figure(2);
imshow(J);

saveas(2, '../../img/2 - Rectification/PalazzoTe_rectified.png');
