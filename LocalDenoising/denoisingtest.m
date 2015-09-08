clear;close all;clc;

addpath(genpath('../'));

%% load image
image = imread('barb.bmp');
image = im2double(image);

%% Add noise
sigmaabs = 20;
sigmasquare = sigmaabs^2/(255^2);
sigma = sqrt(sigmasquare);
noisedimage = imnoise(image,'gaussian',0,sigmasquare);

%% Denoise
denoisedimage = local_Denoise(image,noisedimage,sigma,1);
denoisedimage2 = local_Denoise(image,noisedimage,sigma,2);
% denoisedimage3 = BSdenoise(noisedimage);

%% Calculate MSE and PSNR
MSE = sum(sum(((denoisedimage-image)*255).^2))/(512^2);
PSNR = 10*log10(255^2/MSE);

MSE2 = sum(sum(((denoisedimage2-image)*255).^2))/(512^2);
PSNR2 = 10*log10(255^2/MSE2);

% MSE3 = sum(sum(((denoisedimage3-image)*255).^2))/(512^2);
% PSNR3 = 10*log10(255^2/MSE3);

%% Show image
% figure;
% imshow(image);  % imshow [0 1] grey
% figure;
% imshow(noisedimage);  
% figure;
% imshow(denoisedimage);
% % print('-r600','-dpng','Mandrill-Denoise1');
% figure;
% imshow(denoisedimage2);
