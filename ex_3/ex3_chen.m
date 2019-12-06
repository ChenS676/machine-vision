clear all, close all, clc
L = load('pixellist_postit2g.mat')
for i = 1:1:length(L.pixellist)
    Li = L.pixellist(i).list;
    plot(Li(:,1),Li(:,2), '*b');
    
    
    sum_x = sum(Li(:,1));
    sum_y = sum(Li(:,2));
    sum_x2 = sum(Li(:,1).^2)
    sum_x2 = sum(Li(:,2).^2)
    sum_xy = sum(Li(:,2).*Li(:,1))
    
    
    N = length(Li)
    alpha = sum_x2 - sum_x^2/N
    beta = sum_xy - sum_z*sum_y/N
    gamma = sum_y2 - sum_y.^2/N
    
    
    M = [alpha, beta; beta, gamma];
    [n, ~] = eig(M);
    n = n(:,1);
    c = -sum(n'*Li')/N;
    
    lines(i, :) = [n', c'];
end

x = 1:0.1:size(I,2);
for i = 1: 1: size(lines,1)
    y = (-lines(i,1)*x - lines(i, 3))/lines(i, 2);
    plot(x, y, '-k', 'LineWidth', l);
end



%% Exercise 2
clear all, close all, clc
L = load('pixellist_postit2g.mat')
for i = 1:1:length(L.pixellist)
    Li = L.pixellist(i).list;
    plot(Li(:,1),Li(:,2), '*b');
    
    
    sum_x = sum(Li(:,1));
    sum_y = sum(Li(:,2));
    sum_x2 = sum(Li(:,1).^2)
    sum_x2 = sum(Li(:,2).^2)
    sum_xy = sum(Li(:,2).*Li(:,1))
    
    
    N = length(Li)
    alpha = sum_x2 - sum_x^2/N
    beta = sum_xy - sum_z*sum_y/N
    gamma = sum_y2 - sum_y.^2/N
    
    
    M = [alpha, beta; beta, gamma];
    [n, ~] = eig(M);
    n = n(:,1);
    c = -sum(n'*Li')/N;
    
    lines(i, :) = [n', c'];
    
    d = lines(i, 1)*Li(:,1) + lines(i, 2) *Li(:,2)+lines(i,3);
    Li_projection = Li- [lines(i,1)*d, lines(i,2)*d];
    plot(Li_projection(:, 1), Li_projection(:,2), 'r*', 'Markersize')
    
    tau = Li*lines(i, 1:2);
    
end

x = 1:0.1:size(I,2);
for i = 1: 1: size(lines,1)
    y = (-lines(i,1)*x - lines(i, 3))/lines(i, 2);
    plot(x, y, '-k', 'LineWidth', l);
end