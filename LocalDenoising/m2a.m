function [ image ] = m2a( scale1m, scale2m, scale3m, LL )
%m2a Wavelet transform coefficients matrix to array
%   [ image ] = m2a( scale1m, scale2m, scale3m )
%   scale1m, scale2m, scale3m are three dimensional matrix 
%   "image" is the inverse transform result


scale1 = zeros(3,65536);
scale2 = zeros(3,16384);
scale3 = zeros(3,4096);

for i = 1:3
    for j = 1:256
        scale1(i,(j-1)*256+1:j*256) = scale1m(j,:,i);
    end
end

for i = 1:3
    for j = 1:128
        scale2(i,(j-1)*128+1:j*128) = scale2m(j,:,i);
    end
end

for i = 1:3
    for j = 1:64
        scale3(i,(j-1)*64+1:j*64) = scale3m(j,:,i);
    end
end

C = [LL scale3(1,:) scale3(2,:) scale3(3,:) scale2(1,:) scale2(2,:) scale2(3,:) scale1(1,:) scale1(2,:) scale1(3,:)];
S = [64 64;64 64; 128 128;256 256;512 512];

image = waverec2(C,S,'haar');
    
end

