function [n, c] = fit_line_ls(pixellist)
    N = length(pixellist);
    sum_x = sum(pixellist(:, 1));
    sum_x2 = sum(pixellist(:, 1).^2);
    sum_y = sum(pixellist(:, 2));
    sum_y2 = sum(pixellist(:, 2).^2);
    sum_xy = sum(pixellist(:, 1).*pixellist(:, 2));
    
    alpha = sum_x2 - sum_x^2 / N;
    beta = sum_xy - sum_x * sum_y / N;
    gamma = sum_y2 - sum_y^2 / N;
    
    M = [alpha, beta; beta, gamma];
    [n, ~] = eig(M);
    n = n (:, 1);
    c = - sum(n' * pixellist') / N;
end