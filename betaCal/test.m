clear;clc;

addpath('../BandDenoising');

load('rtable.mat');
image = imread('lena.bmp');
image_double = im2double(image);

[scale1, scale2, scale3] = dwt2_3level(image_double);

HH1 = scale1(1,:);
% mux = sum(HH1)/65536;
% sigmax2 = sum((HH1-mux).^2)/65536;
% E2 = (sum(abs(HH1-mux))/65536)^2;
% rhat = sigmax2/E2;
betahat = Calbeta(HH1,r,beta);