%% assignment 02
clear all;


%% Exercise 1
close all; clc;

%% Load image
I = imread('postit2g.png');
I_d = im2double(I);

%% Sobel
sobel_fm_u = (1/8) * [1, 0, -1;
                   2, 0, -2;
                   1, 0, -1];
sobel_fm_v = sobel_fm_u';
prewitt_fm_u = (1/6) * [1, 0, -1;
                     1, 0, -1;
                     1, 0, -1];
prewitt_fm_v = prewitt_fm_u';

g_u_sobel = conv2(I_d, sobel_fm_u);
g_v_sobel = conv2(I_d, sobel_fm_v);
figure, 
imshow(g_u_sobel)
figure, 
imshow(g_v_sobel)
% ???????????????????
%% 
% flip ??
g_u_sobel_alt = imfilter(I_d, flip(sobel_fm_u, 2));
g_v_sobel_alt = imfilter(I_d, flip(sobel_fm_v, 1));

figure;
subplot(3,1,1);
imshow(I_d);
subplot(3,1,2);
imshow(g_u_sobel, [], 'Colormap', cool);
colorbar;
subplot(3,1,3);
imshow(g_u_sobel_alt, [], 'Colormap', cool);
colorbar;

figure;
subplot(3,1,1);
imshow(I_d);
subplot(3,1,2);
imshow(g_v_sobel, [], 'Colormap', cool);
colorbar;
subplot(3,1,3);
imshow(g_v_sobel_alt, [], 'Colormap', cool);
colorbar;

%% Prewitt
g_u_prewitt = conv2(I_d, prewitt_fm_u);
g_v_prewitt = conv2(I_d, prewitt_fm_v);
 
%% Gradient length and angle
g_sobel_length = sqrt(g_u_sobel.^2 + g_v_sobel.^2);
g_sobel_angle = atan2(g_v_sobel, g_u_sobel);
g_sobel_angle_scaled = (g_sobel_angle / (2 * pi)) + 0.5;
g_prewitt_length = sqrt(g_u_prewitt.^2 + g_v_prewitt.^2);
g_prewitt_angle = atan2(g_v_prewitt, g_u_prewitt);
g_prewitt_angle_scaled = (g_prewitt_angle / (2 * pi)) + 0.5;

figure;
subplot(2,1,1)
imshow(g_sobel_length, [], 'Colormap', jet);
subplot(2,1,2)
imshow(g_prewitt_length, [], 'Colormap', jet);

figure;
subplot(3,1,1);
imshow(I_d);
subplot(3,1,2);
imshow(g_sobel_length, []);
colorbar;
subplot(3,1,3);
imshow(g_prewitt_angle * 180 / pi, [-180, 180], 'Colormap', hsv);
colorbar;
 
%% Exercise 2
%close all; clc;

%% Laplacian and Sobel
laplacian_fm = fspecial('laplacian');
I_laplacian = conv2(I_d, laplacian_fm);
I_laplacian_bw = im2bw(I_laplacian, 0.15);
I_sobel_bw = im2bw(g_sobel_length, 0.05);

figure;
subplot(2,2,1);
imshow(I_laplacian, []);
colorbar;
subplot(2,2,2);
imshow(g_sobel_length, []);
colorbar;
subplot(2,2,3);
imshow(I_laplacian_bw);
subplot(2,2,4);
imshow(I_sobel_bw);
%??sobel???laplacian?robust


%% Canny and Laplacian of Gaussian
% %      The Canny method finds edges by looking for local maxima of the
% %      gradient of I. The gradient is calculated using the derivative of a
% %      Gaussian filter. The method uses two thresholds, to detect strong
% %      and weak edges, and includes the weak edges in the output only if
% %      they are connected to strong edges. This method is therefore less
% %      likely than the others to be "fooled" by noise, and more likely to
% %      detect true weak edges.
I_canny = edge(I_d, 'canny', [0.03, 0.15], 2);


%   BW = EDGE(I,'log',THRESH) specifies the sensitivity threshold for the
%   Laplacian of Gaussian method. EDGE ignores all edges that are not
%   stronger than THRESH. If you do not specify THRESH, or if THRESH is
%   empty ([]), EDGE chooses the value automatically.
%
%   BW = EDGE(I,'log',THRESH,SIGMA) specifies the Laplacian of Gaussian
%   method, using SIGMA as the standard deviation of the LoG filter. The
%   default SIGMA is 2; the size of the filter is N-by-N, where
%   N=CEIL(SIGMA*3)*2+1.
I_log = edge(I_d, 'log', 0.0008, 3);

figure;
subplot(2,1,1);
imshow(I_canny);
subplot(2,1,2);
imshow(I_log);
 
%% Exercise 3
close all; clc;
n_lines = 10;

hough_data = robust_hough(I_canny);
hough_lines = robust_hough_lines(hough_data, n_lines, I_canny);
robust_hough_plot_lines(I_d, hough_lines);
 
[~, index]= sort(hough_data.peaks(:, 3));
peaks = hough_data.peaks(index, :);
peaks = peaks(end - n_lines:end, :);

figure;
imagesc(hough_data.accumulator);
colormap hot;
hold on;
plot(peaks(:,2),peaks(:,1),'wo');
hold off;