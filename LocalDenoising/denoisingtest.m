clear;close all;clc;

addpath('../.');

%% load image
image = imread('lena.bmp');
image = im2double(image);

%% Add noise
sigmaabs = 20;
sigmasquare = sigmaabs^2/(255^2);
sigma = sqrt(sigmasquare);
noisedimage = imnoise(image,'gaussian',0,sigmasquare);

%% Denoise
denoisedimage = local_Denoise(image,noisedimage,sigma,1);
denoisedimage2 = local_Denoise(image,noisedimage,sigma,2);

%% Calculate MSE and PSNR
MSE = sum(sum(((denoisedimage-image)*255).^2))/(512^2);
PSNR = 10*log10(255^2/MSE);

MSE2 = sum(sum(((denoisedimage2-image)*255).^2))/(512^2);
PSNR2 = 10*log10(255^2/MSE2);

%% Show image
figure;
imshow(image);  % imshow [0 1] grey
% print('-r600','-dpng','Mandrill-Original');
figure;
imshow(noisedimage);  
% print('-r600','-dpng','Mandrill-Addnoise');
figure;
imshow(denoisedimage);
% print('-r600','-dpng','Mandrill-Denoise1');
figure;
imshow(denoisedimage2);