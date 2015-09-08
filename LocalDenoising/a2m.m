function [ scale1m, scale2m, scale3m, LL ] = a2m( image )
%a2m Wavelet transform coefficients array to matrix
%   [ scale1m, scale2m, scale3m ] = a2m( image )
%   scale1m, scale2m, scale3m are three dimensional matrix 
%   "image" is the input image to be transformed. 

[ scale1, scale2, scale3 ] = dwt2_3level( image );
scale1m = zeros(256,256,3);
scale2m = zeros(128,128,3);
scale3m = zeros(64,64,3);

LL = scale3(1,:);
for i = 1:3
    for j = 1:256
        scale1m(:,j,i) = scale1(i,(j-1)*256+1:j*256);
    end
end

for i = 1:3
    for j = 1:128
        scale2m(:,j,i) = scale2(i,(j-1)*128+1:j*128);
    end
end

for i = 1:3
    for j = 1:64
        scale3m(:,j,i) = scale3(i+1,(j-1)*64+1:j*64);
    end
end

    
end

