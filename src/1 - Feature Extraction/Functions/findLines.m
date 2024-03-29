function lines = findLines(I)
% HOUGH TRANSFORM
    
    % Canny Algorithm
    edges = edge(I, 'canny', [0.2 0.25]);
    
    % Hough Transform space
    [H,T,R] = hough(edges);
    
    % Peaks in Hough Transform matrix
    P = houghpeaks(H, 100, 'threshold', 0.3*max(H(:)));
    
    % Extract line segments
    lines = houghlines(edges, T, R, P, 'FillGap', 200, 'MinLength', 500);
end
