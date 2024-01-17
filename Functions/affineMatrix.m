function [H, imLinf] = buildHaff(lines)
    V = nan(2, length(lines));
    for i = 1:length(lines)
        A = lines{i}(:,1:2);
        B = -lines{i}(:,3);
        V(:,i) = A\B;
    end
    imLinf = fitLine(V);
    H = [eye(2), zeros(2,1); imLinf(:).'];
end