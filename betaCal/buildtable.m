clear;clc;

beta = 0.001:0.001:10;
r = gamma(1./beta).*gamma(3./beta)./(gamma(2./beta).^2);
save('rtable.mat','beta','r');