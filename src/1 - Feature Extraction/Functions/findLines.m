function lines = findLines(I)
% HOUGH TRANSFORM
    
    % Canny Algorithm
    edges = edge(I, 'canny', [0.2 0.25]);
    
    % Hough Transform space
    [H,T,R] = hough(edges, 'RhoResolution', 0.5, 'Theta', -90:0.5:89.5);
    
    % Peaks in Hough Transform matrix
    P = houghpeaks(H, 50, 'threshold', ceil(0.3*max(H(:))), 'NHoodSize', [97 37]);
    
    % Extract line segments
    lines = houghlines(edges, T, R, P, 'FillGap', 200, 'MinLength', 500);
end
