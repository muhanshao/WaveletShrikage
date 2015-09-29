clear;close all;clc;

addpath(genpath('../'));

load('rtable.mat');
image = imread('barb.bmp');
image = im2double(image);

[ sigmax1m, sigmax2m, sigmax3m, beta1m, beta2m, beta3m ] = local_cal( image,4,r,beta );


% Level 1
beta11 = beta1m(:,:,1);
figure;
subplot(331);
histogram(beta11);
title('LH');
xlabel('\beta');
ylabel('Level 1');

beta12 = beta1m(:,:,2);
subplot(332);
histogram(beta12);
title('HL');
xlabel('\beta');

beta13 = beta1m(:,:,3);
subplot(333);
histogram(beta13);
title('HH');
xlabel('\beta');

% Level 2
beta21 = beta2m(:,:,1);
subplot(334);
histogram(beta21);
title('LH');
xlabel('\beta');
ylabel('Level 2');

beta22 = beta2m(:,:,2);
subplot(335);
histogram(beta22);
title('HL');
xlabel('\beta');

beta23 = beta2m(:,:,3);
subplot(336);
histogram(beta23);
title('HH');
xlabel('\beta');

% Level 3
beta31 = beta3m(:,:,1);
subplot(337);
histogram(beta31);
title('LH');
xlabel('\beta');
ylabel('Level 3');

beta32 = beta3m(:,:,2);
subplot(338);
histogram(beta32);
title('HL');
xlabel('\beta');

beta33 = beta3m(:,:,3);
subplot(339);
histogram(beta33);
title('HH');
xlabel('\beta');


% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'points');
% set(gcf, 'PaperPosition', [0 0 1000 1000]);
% print('-r300','-dpng','Histogram_barb');