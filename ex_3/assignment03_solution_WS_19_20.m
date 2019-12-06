%% assignment 03
clear all;
close all; 
clc;

%% Exercise 1
I = imread('postit2g.png');
L = load('pixellist_postit2g.mat');

figure;
imagesc(I);
colormap gray;
hold on;

for i = 1:1:length(L.pixellist)
    Li = L.pixellist(i).list;
    plot(Li(:,1), Li(:,2), '*b', 'LineWidth', 3);
    
    sum_x = sum(Li(:, 1));
    sum_x2 = sum(Li(:, 1).^2);
    sum_y = sum(Li(:, 2));
    sum_y2 = sum(Li(:, 2).^2);
    sum_xy = sum(Li(:, 1).*Li(:, 2));
    
    alpha = sum_x2 - sum_x^2 / size(Li, 1);
    beta = sum_xy - sum_x * sum_y / size(Li, 1);
    gamma = sum_y2 - sum_y^2 / size(Li, 1);
    
    M = [alpha, beta; beta, gamma];
    [n, ~] = eig(M);
    
    % minimum corresponds to the first eigenvector
    n = n(:, 1);
    c = - sum(n' * Li') / size(Li, 1);
    
    lines(i, :) = [n', c];
    
    x = 1:0.1:size(I,2);
    y = (-lines(i, 1) * x - lines(i, 3)) / lines(i, 2);
    plot(x, y, '-k', 'LineWidth', 1);
end


%% Exercise 2
for i = 1:1:length(L.pixellist)
    Li = L.pixellist(i).list;
    d = lines(i, 1) *  Li(:, 1) + lines(i, 2) *  Li(:, 2) + lines(i, 3);
    Li_proj = Li - [lines(i, 1) * d, lines(i, 2) * d];    
    plot(Li_proj(:,1), Li_proj(:,2), '*r', 'LineWidth', 3);
    
    % minimum and maximum in x direction
    [~, min_i] = min(Li_proj(:,1));
    [~, max_i] = max(Li_proj(:,1));  
    
    % minimum and maximum orthogonal to the normal vector
    tau = Li * [-lines(i, 2); lines(i, 1)];
    [~, tau_min_i] = min(tau);
    [~, tau_max_i] = max(tau);
    
    plot([Li_proj(min_i,1), Li_proj(max_i,1)], [Li_proj(min_i,2), Li_proj(max_i,2)], '-og', 'LineWidth', 3);
    plot([Li_proj(tau_min_i,1), Li_proj(tau_max_i,1)], [Li_proj(tau_min_i,2), Li_proj(tau_max_i,2)], '-ob', 'LineWidth', 1);
end

hold off;

%% Exercise 3
img = figure;
imagesc(I);
colormap gray;
hold on;

Li = [L.pixellist(1).list; 
      L.pixellist(2).list];
plot(Li(:,1), Li(:,2), '*b', 'LineWidth', 2);

% pick 2 points "ramdomly"
perm = [round(size(Li, 1) * 0.075), round(size(Li,1) * 0.6)];
plot(Li(perm,1), Li(perm,2), 'xc', 'MarkerSize', 10);
% threshold: maximum tolerable error until convergence is reached
error = 10^(-2);
% parameter for the animation
duration = 1;

%% LS
[n, c] = fit_line_ls(Li);

ls_start_end = find_start_end(n, c, Li);
plot(ls_start_end(:,1), ls_start_end(:,2), '-og', 'LineWidth', 2);
lgd = legend('Measured Points','Random Points', 'Least Squares');
lgd.FontSize = 16;
pause(duration);

%% M-Estimators
x = -30:0.1:30;
y = x.^2;
dy = 2*x;

k = 3;
Huber(x >= -k & x <= k) = x(x >= -k & x <= k).^2;
Huber(x < -k | x > k) = 2 * k * abs(x(x < -k | x > k)) - k^2;
dHuber(x >= -k & x <= k) = 2 * x(x >= -k & x <= k);
dHuber(x < -k) = -2 * k; 
dHuber(x > k) = 2 * k;

kappa = 3;
Cauchy = kappa^2 * log(1 + x.^2/kappa^2);
dCauchy = 2 * x ./ (1 + x.^2/kappa^2);

m_estimators = figure;
plot(x, y, '-b', 'LineWidth', 3);
hold on;
axis([-5, 5, -10, 30]);
plot(x, Huber, '-g', 'LineWidth', 3);
plot(x, Cauchy, '-r', 'LineWidth', 3);
plot(x, dy, '--b', 'LineWidth', 3);
plot(x, dHuber, '--g', 'LineWidth', 3);
plot(x, dCauchy, '--r', 'LineWidth', 3);
lgd2 = legend('Least Squares','Huber','Cauchy');
lgd2.FontSize = 16;

%% M-Estimators - Huber
subset = [Li(perm(1),:);
          Li(perm(2),:)];

% calculate line parameters for 2 points
dxy = Li(perm(1),:) - Li(perm(2),:);
n = [-dxy(2); dxy(1)] ./sqrt(dxy*dxy');
c = -n'*Li(perm(1), :)';

d = (n(1) * Li(:, 1) + n(2) * Li(:, 2) + c);
dn = [1; 1];
dc = 1;
    
huber_start_end = find_start_end(n, c, Li);
figure(img);
huber_plot(1) = plot(huber_start_end(:,1), huber_start_end(:,2), '--oy', 'LineWidth', 1);
pause(duration);

i = 2;
while (dn(1) > error || dn(2) > error || dc > error)
    % calculate weights
    weights(d >= -k & d <= k) = 1;
    weights(d < -k) = -2 * k./(2 * d(d < -k));
    weights(d > k) = 2 * k./(2 * d(d > k));
    
    n_before = n;
    c_before = c;
    [n, c] = fit_line_wls(Li, weights');
    dn = n - n_before;
    dc = c - c_before;
    
    huber_start_end = find_start_end(n, c, Li);
    figure(img);
    huber_plot(i) = plot(huber_start_end(:,1), huber_start_end(:,2), '--oy', 'LineWidth', 1);
    pause(duration);
    i = i + 1;
    
    % calculate d to get a new distribution of the weights in the next iteration
    d = (n(1) * Li(:, 1) + n(2) * Li(:, 2) + c);
end

huber_start_end = find_start_end(n, c, Li);
figure(img);
plot(huber_start_end(:,1), huber_start_end(:,2), '-oy', 'LineWidth', 2);
delete(huber_plot);
lgd = legend('Measured Points','Random Points','Least Squares','M-Estimator + Huber');
lgd.FontSize = 16;
pause(duration);

%% M-Estimators - Cauchy
subset = [Li(perm(1),:);
          Li(perm(2),:)];
      
% calculate line parameters for 2 points
dxy = Li(perm(1),:) - Li(perm(2),:);
n = [-dxy(2); dxy(1)] ./sqrt(dxy*dxy');
c = -n'*Li(perm(1), :)';
    
d = (n(1) * Li(:, 1) + n(2) * Li(:, 2) + c);
dn = [1; 1];
dc = 1;

cauchy_start_end = find_start_end(n, c, Li);
figure(img);
cauchy_plot(1) = plot(cauchy_start_end(:,1), cauchy_start_end(:,2), '--or', 'LineWidth', 1);
pause(duration);

i = 2;
while (dn(1) > error || dn(2) > error || dc > error)
    % calculate weights
    weights = 1./(1 + (d/kappa).^2);
    
    n_before = n;
    c_before = c;
    [n, c] = fit_line_wls(Li, weights);
    dn = n - n_before;
    dc = c - c_before;
    
    cauchy_start_end = find_start_end(n, c, Li);
    figure(img);
    cauchy_plot(i) = plot(cauchy_start_end(:,1), cauchy_start_end(:,2), '--or', 'LineWidth', 1);
    pause(duration);
    i = i + 1;
        
    % calculate d to get a new distribution of the weights in the next iteration
    d = (n(1) * Li(:, 1) + n(2) * Li(:, 2) + c);
end

cauchy_start_end = find_start_end(n, c, Li);
figure(img);
plot(cauchy_start_end(:,1), cauchy_start_end(:,2), '-or', 'LineWidth', 2);
delete(cauchy_plot);
lgd = legend('Measured Points','Random Points','Least Squares','M-Estimator + Huber','M-Estimator + Cauchy');
lgd.FontSize = 16;
pause(duration);

%% LTS
% acceptance rate
p = 0.8;
p_points = round(size(Li, 1) * p);

subset = [Li(perm(1),:);
          Li(perm(2),:)];

% calculate line parameters for 2 points
dxy = Li(perm(1),:) - Li(perm(2),:);
n = [-dxy(2); dxy(1)] ./sqrt(dxy*dxy');
c = -n'*Li(perm(1), :)';

d = (n(1) * Li(:, 1) + n(2) * Li(:, 2) + c);
[~, d_i]= sort(d);
subset = Li(d_i(1:p_points), :);
dn = [1; 1];
dc = 1;

lts_start_end = find_start_end(n, c, Li);
figure(img);
lts_plot(1) = plot(lts_start_end(:,1), lts_start_end(:,2), '--om', 'LineWidth', 1);
pause(duration);

i = 2;
while (dn(1) > error || dn(2) > error || dc > error)
    n_before = n;
    c_before = c; 
    [n, c] = fit_line_ls(subset);
    dn = n - n_before;
    dc = c - c_before;

    lts_start_end = find_start_end(n, c, Li);
    figure(img);
    lts_plot(i) = plot(lts_start_end(:,1), lts_start_end(:,2), '--om', 'LineWidth', 1);
    pause(duration);
    i = i + 1;
    
    % calculate d and sort points to determine a new subset
    d = abs((n(1) * Li(:, 1) + n(2) * Li(:, 2) + c));
    [~, d_i]= sort(d);
    subset = Li(d_i(1:p_points), :);
end

lts_start_end = find_start_end(n, c, Li);
figure(img);
plot(lts_start_end(:,1), lts_start_end(:,2), '-om', 'LineWidth', 2);
delete(lts_plot);
lgd = legend('Measured Points','Random Points','Least Squares','M-Estimator + Huber','M-Estimator + Cauchy','Least Trimmed Squares');
lgd.FontSize = 16;
pause(duration);

%% RANSAC
% threshold: tolerable distance to the fitted line
thresh = 2;
% number of trials
i_max = 20;

for i = 1:1:i_max
    % pick 2 points randomly
    perm = randperm(size(Li, 1), size(Li, 2));
    
    % calculate line parameters for 2 points
    dxy = Li(perm(1),:) - Li(perm(2),:);
    n = [-dxy(2); dxy(1)] ./sqrt(dxy*dxy');
    c = -n'*Li(perm(1), :)';
    
    ransac_start_end = find_start_end(n, c, Li);
    figure(img);
    ransac_plot(i)= plot(ransac_start_end(:,1), ransac_start_end(:,2), '--ok', 'LineWidth', 1);
    pause(duration);
    
    % check how many points are outside the tolerance band 
    d = sum(abs((n(1) * Li(:, 1) + n(2) * Li(:, 2) + c)) > thresh);
    trials(i, :) = [n', c, d];
end

[~, trials_i] = sort(trials(:, 4));
ransac_start_end = find_start_end(trials(trials_i(1), 1:2)', trials(trials_i(1), 3), Li);
figure(img);
plot(ransac_start_end(:,1), ransac_start_end(:,2), '-ok', 'LineWidth', 2);
delete(ransac_plot);

lgd = legend('Measured Points','Random Points','Least Squares','M-Estimator + Huber','M-Estimator + Cauchy','Least Trimmed Squares','RANSAC');
lgd.FontSize = 16;
pause(duration);


