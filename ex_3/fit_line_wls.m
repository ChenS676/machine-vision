function [n, c] = fit_line_wls(pixellist, weights)
    W = sum(weights);
    sum_x = sum(weights.*pixellist(:, 1));
    sum_x2 = sum(weights.*pixellist(:, 1).^2);
    sum_y = sum(weights.*pixellist(:, 2));
    sum_y2 = sum(weights.*pixellist(:, 2).^2);
    sum_xy = sum(weights.*pixellist(:, 1).*pixellist(:, 2));
    
    alpha = sum_x2 - sum_x^2 / W;
    beta = sum_xy - sum_x * sum_y / W;
    gamma = sum_y2 - sum_y^2 /W;
    
    M = [alpha, beta; beta, gamma];
    [n, ~] = eig(M);
    n = n (:, 1);
    c = - sum(n' * (weights.*pixellist)') / W;
end