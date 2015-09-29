clear;close all;clc;

addpath(genpath('../'));

%% load image
image = imread('lena.bmp');
image = im2double(image);

%% Add noise
sigmaabs = 100;
sigmasquare = sigmaabs^2/(255^2);
sigma = sqrt(sigmasquare);
% noisedimage = imnoise(image,'gaussian',0,sigmasquare);

testnum = 5;
MSE = zeros(testnum,16);

%% Denoise
for i = 1:testnum
    noisedimage = imnoise(image,'gaussian',0,sigmasquare);
    % local denoise
    % soft-thresholding
    % real sigmax, the original function
    denoisedimage1 = local_Denoise(image,noisedimage,sigma,1,'s'); 
    % real sigmax, the modified function
    denoisedimage2 = local_Denoise(image,noisedimage,sigma,2,'s');  
    % estimated sigmax, the original function
    denoisedimage3 = local_Denoise(image,noisedimage,sigma,3,'s'); 
    % estimated sigmax, the modified function
    denoisedimage4 = local_Denoise(image,noisedimage,sigma,4,'s'); 
    % hard thresholding
    denoisedimage5 = local_Denoise(image,noisedimage,sigma,1,'h');
    denoisedimage6 = local_Denoise(image,noisedimage,sigma,2,'h');
    denoisedimage7 = local_Denoise(image,noisedimage,sigma,3,'h');
    denoisedimage8 = local_Denoise(image,noisedimage,sigma,4,'h');
  
    MSE(i,1) = sum(sum(((denoisedimage1-image)*255).^2))/(512^2);
    MSE(i,2) = sum(sum(((denoisedimage2-image)*255).^2))/(512^2);
    MSE(i,3) = sum(sum(((denoisedimage3-image)*255).^2))/(512^2);
    MSE(i,4) = sum(sum(((denoisedimage4-image)*255).^2))/(512^2);
    MSE(i,5) = sum(sum(((denoisedimage5-image)*255).^2))/(512^2);
    MSE(i,6) = sum(sum(((denoisedimage6-image)*255).^2))/(512^2);
    MSE(i,7) = sum(sum(((denoisedimage7-image)*255).^2))/(512^2);
    MSE(i,8) = sum(sum(((denoisedimage8-image)*255).^2))/(512^2);
    
    % band denoise
    % soft-thresholding
    % real sigmax, the original function
    denoisedimage11 = band_Denoise(image,noisedimage,sigma,1,'s'); 
    % real sigmax, the modified function
    denoisedimage22 = band_Denoise(image,noisedimage,sigma,2,'s');  
    % estimated sigmax, the original function
    denoisedimage33 = band_Denoise(image,noisedimage,sigma,3,'s'); 
    % estimated sigmax, the modified function
    denoisedimage44 = band_Denoise(image,noisedimage,sigma,4,'s'); 
    % hard thresholding
    denoisedimage55 = band_Denoise(image,noisedimage,sigma,1,'h');
    denoisedimage66 = band_Denoise(image,noisedimage,sigma,2,'h');
    denoisedimage77 = band_Denoise(image,noisedimage,sigma,3,'h');
    denoisedimage88 = band_Denoise(image,noisedimage,sigma,4,'h');
  
    MSE(i,9) = sum(sum(((denoisedimage11-image)*255).^2))/(512^2);
    MSE(i,10) = sum(sum(((denoisedimage22-image)*255).^2))/(512^2);
    MSE(i,11) = sum(sum(((denoisedimage33-image)*255).^2))/(512^2);
    MSE(i,12) = sum(sum(((denoisedimage44-image)*255).^2))/(512^2);
    MSE(i,13) = sum(sum(((denoisedimage55-image)*255).^2))/(512^2);
    MSE(i,14) = sum(sum(((denoisedimage66-image)*255).^2))/(512^2);
    MSE(i,15) = sum(sum(((denoisedimage77-image)*255).^2))/(512^2);
    MSE(i,16) = sum(sum(((denoisedimage88-image)*255).^2))/(512^2);
   
end

PSNR = 10*log10(255^2./MSE);


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
