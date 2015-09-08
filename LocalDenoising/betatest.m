clear;close all;clc;

addpath(genpath('../'));

load('rtable.mat');
image = imread('lena.bmp');
image = im2double(image);

[ sigmax1m, sigmax2m, sigmax3m, beta1m, beta2m, beta3m ] = local_cal( image,4,r,beta );

beta11 = beta1m(:,:,1);
figure;
subplot(131);
histogram(beta11);
beta12 = beta1m(:,:,2);
subplot(132);
histogram(beta12);
beta13 = beta1m(:,:,3);
subplot(133);
histogram(beta13);

beta21 = beta2m(:,:,1);
figure;
subplot(131);
histogram(beta21);
beta22 = beta2m(:,:,2);
subplot(132);
histogram(beta22);
beta23 = beta2m(:,:,3);
subplot(133);
histogram(beta23);

beta31 = beta3m(:,:,1);
figure;
subplot(131);
histogram(beta31);
beta32 = beta3m(:,:,2);
subplot(132);
histogram(beta32);
beta33 = beta3m(:,:,3);
subplot(133);
histogram(beta33);

