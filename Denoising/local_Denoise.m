function [ denoisedimage ] = local_Denoise( image, noisedimage, sigma, method, thrtype )
%% local_Denoise Denoise a image using wavelet thresholding method and calculate sigmax locally.
%   [denoisedimage] = local_Denoise(image, noisedimage, sigma, method, thrtype);
%   image: clean image
%   noisedimage: noised image
%   sigma: noise standard deviation
%   method = 1: calculate real sigmax using clean image, thresholding by original formula
%   method = 2: calculate real sigmax using clean image, thresholding by modified formula
%   method = 3: estimate sigmax using noised image, thrsholding by original formula
%   method = 4: estimate sigmax using noised image, thrsholding by modified formula
%   thrtype: thresholding type. 'h': hard-thresholding. 's': soft-thresholding.

[scale1m, scale2m, scale3m, LL] = a2m(noisedimage);

descale1m = zeros(256,256,3);
descale2m = zeros(128,128,3);
descale3m = zeros(64,64,3);

% cons = 1.25;

switch method
    case {1}   % clean image, original formula
        [sigmax1m, sigmax2m, sigmax3m, flag1m, flag2m, flag3m] = local_sigmax(image,4,sigma,1);
        for i = 1:3
            for j = 1:256
                for k = 1:256
                    tempT = Tbo(sigma,sigmax1m(j,k,i),flag1m(j,k,i));
                    descale1m(j,k,i) = wthresh(scale1m(j,k,i),thrtype,tempT);
                end
            end
        end
        for i = 1:3
            for j = 1:128
                for k = 1:128
                    tempT = Tbo(sigma,sigmax2m(j,k,i),flag2m(j,k,i));
                    descale2m(j,k,i) = wthresh(scale2m(j,k,i),thrtype,tempT);
                end
            end
        end
        for i = 1:3
            for j = 1:64
                for k = 1:64
                    tempT = Tbo(sigma,sigmax3m(j,k,i),flag3m(j,k,i));
                    descale3m(j,k,i) = wthresh(scale3m(j,k,i),thrtype,tempT);
                end
            end
        end
      
    case {2}  % clean image, modified function
        [sigmax1m, sigmax2m, sigmax3m, flag1m, flag2m, flag3m] = local_sigmax(image,4,sigma,1);
        for i = 1:3
            for j = 1:256
                for k = 1:256
                    tempT = Tb051(sigma,sigmax1m(j,k,i),flag1m(j,k,i));
                    descale1m(j,k,i) = wthresh(scale1m(j,k,i),thrtype,tempT);
                end
            end
        end
        for i = 1:3
            for j = 1:128
                for k = 1:128
                    tempT = Tbo(sigma,sigmax2m(j,k,i),flag2m(j,k,i));
                    descale2m(j,k,i) = wthresh(scale2m(j,k,i),thrtype,tempT);
                end
            end
        end
        for i = 1:3
            for j = 1:64
                for k = 1:64
                    tempT = Tbo(sigma,sigmax3m(j,k,i),flag3m(j,k,i));
                    descale3m(j,k,i) = wthresh(scale3m(j,k,i),thrtype,tempT);
                end
            end
        end
     
    case {3}
        [sigmax1m, sigmax2m, sigmax3m, flag1m, flag2m, flag3m] = local_sigmax(noisedimage,4,sigma,2);
        for i = 1:3
            for j = 1:256
                for k = 1:256
                    tempT = Tbo(sigma,sigmax1m(j,k,i),flag1m(j,k,i));
                    descale1m(j,k,i) = wthresh(scale1m(j,k,i),thrtype,tempT);
                end
            end
        end
        for i = 1:3
            for j = 1:128
                for k = 1:128
                    tempT = Tbo(sigma,sigmax2m(j,k,i),flag2m(j,k,i));
                    descale2m(j,k,i) = wthresh(scale2m(j,k,i),thrtype,tempT);
                end
            end
        end
        for i = 1:3
            for j = 1:64
                for k = 1:64
                    tempT = Tbo(sigma,sigmax3m(j,k,i),flag3m(j,k,i));
                    descale3m(j,k,i) = wthresh(scale3m(j,k,i),thrtype,tempT);
                end
            end
        end
        
    case {4}
        [sigmax1m, sigmax2m, sigmax3m, flag1m, flag2m, flag3m] = local_sigmax(noisedimage,4,sigma,2);
        for i = 1:3
            for j = 1:256
                for k = 1:256
                    tempT = Tb051(sigma,sigmax1m(j,k,i),flag1m(j,k,i));
                    descale1m(j,k,i) = wthresh(scale1m(j,k,i),thrtype,tempT);
                end
            end
        end
        for i = 1:3
            for j = 1:128
                for k = 1:128
                    tempT = Tbo(sigma,sigmax2m(j,k,i),flag2m(j,k,i));
                    descale2m(j,k,i) = wthresh(scale2m(j,k,i),thrtype,tempT);
                end
            end
        end
        for i = 1:3
            for j = 1:64
                for k = 1:64
                    tempT = Tbo(sigma,sigmax3m(j,k,i),flag3m(j,k,i));
                    descale3m(j,k,i) = wthresh(scale3m(j,k,i),thrtype,tempT);
                end
            end
        end
        
    otherwise
        disp('Unknown method!');        
end

denoisedimage = m2a(descale1m, descale2m, descale3m, LL);

end

%% Calculate T
% Original
function [T] = Tbo(sigma,sigmax,flag)

if flag
    T = sigmax;
else
    T = sigma^2/sigmax;
end

end

% Gaussian
function [T] = Tb20(sigma,sigmax,flag)

if flag
    T = sigmax;
else
    x = sigmax/sigma;
    T = sigma^2/sigmax*(0.0944/(x+0.0119)+0.78);
end

end

% Laplacian
function [T] = Tb10(sigma,sigmax,flag)
if flag
    T = sigmax;
else
    x = sigmax/sigma;   
    T = sigma^2/sigmax*(-0.01903*x.^2+0.1642*x+0.8104);
end
end

% Beta = 0.5,method 1(multiple)
function [T] = Tb051(sigma,sigmax,flag)

if flag
    T = sigmax;
else
    x = sigmax/sigma;
% T = sigma^2/sigmax*(0.3228*(x)+0.7691);
    T = sigma^2/sigmax*(-0.06958*x^2+0.6081*x+0.5672);
%     T = sigma^2/sigmax*(-0.0112*x^4+0.1161*x^3-0.4642*x^2+1.09*x+0.4191);
end

end

% Beta = 0.5,method 2(divide)
function [T] = Tb052(sigma,sigmax,flag)

if flag
    T = sigmax;
else
    x = sigmax/sigma;
% T = sigma^2/sigmax/(0.3228*(x)+0.7691);
    T = sigma^2/sigmax/(-0.06958*x^2+0.6081*x+0.5672);
end

end