function [corner_x, corner_y] = findCorners(I, sigma, border_margin)
% HARRIS ALGORITHM

    % Derivative masks
    dx = [-1 0 1; 
          -1 0 1; 
          -1 0 1];
    dy = dx';
      
    % Image derivatives -- applying the filter through convolution
    Ix = conv2(I, dx, 'same');
    Iy = conv2(I, dy, 'same');

    % Gaussian filter
    g = fspecial('gaussian', max(1,fix(3*sigma)+1), sigma);

    % Applying gaussian filter to images
    Ix2 = conv2(Ix.^2, g, 'same');
    Iy2 = conv2(Iy.^2, g, 'same');
    Ixy = conv2(Ix.*Iy, g, 'same');

    % Combining filtered images
    cm = (Ix2.*Iy2 - Ixy.^2)./(Ix2 + Iy2 + eps);

    % Set to 0 near boundaries
    cm(1:border_margin,:) = 0;
    cm(end-border_margin:end,:) = 0;
    cm(:,end-border_margin:end) = 0;
    cm(:,1:border_margin) = 0;

    % Threshold cim
    T = mean(cm(:));
    CIM = cm;
    CIM(cm<T) = 0;

    % Perform nonlocal maximum suppression on thresholded measure
    support = true(11);
    
    % Compute maximum over a square neighbor of size 11 x 11
    maxima = ordfilt2(CIM, sum(support(:)), support);
    
    % Determine locations where max corresponds to cim values
    [corner_x, corner_y] = find((cm==maxima).*(CIM>0));
end
