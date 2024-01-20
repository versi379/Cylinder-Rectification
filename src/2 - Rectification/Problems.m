close all
clear
clc
addpath(genpath('Functions'));
%% READ IMAGE
img = imread("../../img/PalazzoTe.jpg");
img = imrotate(img,270);
img = rgb2gray(img); % converts truecolor image RGB to grayscale image
img = im2double(img); % converts image to double precision
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

% two points on C1 (generatrix lines endpoints) -- homogeneous coordinates
p1 = [x(1), y(1), 1];
p2 = [x(5), y(5), 1];
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

% two points on C2 (generatrix lines endpoints) -- homogeneous coordinates
p3 = [x(1), y(1), 1];
p4 = [x(5), y(5), 1];
%% GENERATRIX LINES L1,L2
L1 = cross(p1,p3);
L2 = cross(p2,p4);
%% IMAGE CONTOUR
syms 'x';
syms 'y';

% C1,C2 equations -- cartesian coordinates
eq1 = a1*x^2 + b1*x*y + c1*y^2 + d1*x + e1*y + f1;
eq2 = a2*x^2 + b2*x*y + c2*y^2 + d2*x + e2*y + f2;
eqns = [eq1 == 0, eq2 == 0];

% plot cross sections and generatrix lines
fimplicit(eqns, 'LineWidth',2, 'Color','b');
plot([p1(1), p3(1)], [p1(2), p3(2)], 'LineWidth',2, 'Color','r');
plot([p2(1), p4(1)], [p2(2), p4(2)], 'LineWidth',2, 'Color','r');

% save image
saveas(1, '../../img/2 - Rectification/PalazzoTe_contour.png');
%% CIRCULAR POINTS IMAGE

% C1,C2 intersections
S = solve(eqns, [x,y]);
s1 = [double(S.x(1));double(S.y(1));1];
s2 = [double(S.x(2));double(S.y(2));1];
s3 = [double(S.x(3));double(S.y(3));1];
s4 = [double(S.x(4));double(S.y(4));1];

% circular points image (complex conjugate pair)
II = s1;
JJ = s2;

% dual conic to circular points image
imDCCP = II*JJ' + JJ*II';
imDCCP = imDCCP./norm(imDCCP);
%% VANISHING LINE

% computation
h = cross(II, JJ);
h = h/h(3);

% plot
x = linspace(1,100000,1000000);
y = ((JJ(2)-II(2))/(JJ(1)-II(1)))*(x-II(1))+II(2);
plot(x,y, 'linewidth',2,'color','g');

% save image
saveas(1, '../../img/2 - Rectification/PalazzoTe_vanishing_line.png');
%% C1,C2 DIAMETERS

% computation
d1 = C1*II;
d1 = d1/d1(3);
d2 = C1*JJ;
d2 = d2/d2(3);
d3 = C2*II;
d3 = d3/d3(3);
d4 = C2*JJ;
d4 = d4/d4(3);
%% C1,C2 CENTRES


% computation
O1 = cross(d1,d2);
O1 = O1./O1(3);
O2 = cross(d3,d4);
O2 = O2./O2(3);

% plot
scatter(O1(1), O1(2), 40,'filled','o','MarkerFaceColor','y');
scatter(O2(1), O2(2), 40,'filled','o','MarkerFaceColor','y');
%% AXIS IMAGE

% computation
a = cross(O1,O2);

% plot
x = linspace(1,100000,1000000);
y=((O2(2)-O1(2))/(O2(1)-O1(1)))*(x-O1(1))+O1(2);
plot(x,y,'linewidth',2,'color','c')

% save image
saveas(1, '../../img/2 - Rectification/PalazzoTe_axis.png');
%% AXIS IMAGE VANISHING POINT
lInf = [0;0;1];
vp = cross(a,lInf);
vp = vp/vp(3);
%% CALIBRATION MATRIX
close all
img = imread("../../img/PalazzoTe.jpg");
img = imrotate(img,270);
figure(1), imshow(img);
hold on;
%% pair of parallel lines to find vertical vanishing point vp1
segment1 = drawline('Color','green');
segment2 = drawline('Color','green');
l1 = segToLine(segment1.Position);
l1 = l1 / l1(3);
l2 = segToLine(segment2.Position);
l2 = l2 / l2(3);
vp1 = cross(l1,l2);
vp1 = vp1 / vp1(3);
%% two pairs of parallel lines to find horizontal vanishing point vp2,vp3
segment3 = drawline('Color','red');
segment4 = drawline('Color','red');
segment5 = drawline('Color','red');
segment6 = drawline('Color','red');
l3 = segToLine(segment3.Position);
l3 = l3 / l3(3);
l4 = segToLine(segment4.Position);
l4 = l4 / l4(3);
l5 = segToLine(segment5.Position);
l5 = l5 / l5(3);
l6 = segToLine(segment6.Position);
l6 = l6 / l6(3);
vp2 = cross(l3,l4);
vp2 = vp2 / vp2(3);
vp3 = cross(l5,l6);
vp3 = vp3 / vp3(3);
%% image of the line at the infinity orthogonal to the vertical vanishing points
% in order to obtain two constraints for the calibration process
inf_line = cross(vp2, vp3);
inf_line = inf_line / inf_line(3);
%% affine rectification matrix
H = [eye(2), zeros(2,1); inf_line(:)'];
%% Image of Absolute Conic
IAC = get_IAC(inf_line, vp1, vp2, vp3, H);
%% intrinsic parameters
alfa = sqrt(IAC(1,1));
u0 = -IAC(1,3)/(alfa^2);
v0 = -IAC(2,3);
fy = sqrt(IAC(3,3) - (alfa^2)*(u0^2) - (v0^2));
fx = fy / alfa;
%% build K using the parametrization
K = [fx 0 u0; 0 fy v0; 0 0 1];
disp(K);
%% CYLINDER AXIS ORIENTATION W.R.T. CAMERA REFERENCE
P = [K [0 0 0]'];
pih = P.'*h;
theta1 = 90 + rad2deg(atan(pih(1)/pih(2)));
theta2 = 90 + rad2deg(atan(pih(1)/pih(3)));
%% RATIO: (CIRCULAR CROSS SECTION) RADIUS / DISTANCE

% compute cross section radius
radiusC1 = pdist([p1(1),p1(2);O1(1),O1(2)], 'euclidean');
radiusC2 = pdist([p3(1),p3(2);O2(1),O2(2)], 'euclidean');

% compute distance between cross sections (centers)
distance = pdist([O1(1),O1(2);O2(1),O2(2)], 'euclidean');

% compute ratio
ratio1 = radiusC1/distance; % 1.0991 (computed result)
ratio2 = radiusC2/distance; % 1.5525 (computed result)
