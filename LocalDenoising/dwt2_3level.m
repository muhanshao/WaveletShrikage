function [ scale1, scale2, scale3 ] = dwt2_3level( image )
%dwt2_3level  3 level 2-D wavelet decomposition
%   [scale1,scale2,scale3] = dwt2_3level(image).
%   Image matrix should be 512*512.
%   Wavelet type is 'haar'.
%   scale1 is a 3*(256*256) matrix.
%   scale2 is 3*(128*128)
%   scale3 is 3*(64*64)
%   Every level contains LH(H), HL(V) and HH(D) coefficients.

C = wavedec2(image,3,'haar');
scale1 = zeros(3,65536);
scale2 = zeros(3,16384);
scale3 = zeros(4,4096);

%% scale3
% scale3 LL
scale3(1,:) = C(1:4096);
% scale3 LH 
scale3(2,:) = C(4097 : 8192);
% scale3 HL
scale3(3,:) = C(8193 : 12288);
% scale3 HH
scale3(4,:) = C(12289 : 16384);

%% scale2
% scale2 LH 
scale2(1,:) = C(16385 : 32768);
% scale2 HL
scale2(2,:) = C(32769 : 49152);
% scale2 HH
scale2(3,:) = C(49153 : 65536);

%% scale1
% scale1 LH 
scale1(1,:) = C(65537 : 131072);
% scale1 HL
scale1(2,:) = C(131073 : 196608);
% scale1 HH
scale1(3,:) = C(196609 : 262144);

end

