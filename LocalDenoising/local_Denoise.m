function [ denoisedimage ] = local_Denoise( image, noisedimage, sigma, method )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[sigmax1m, sigmax2m, sigmax3m] = local_sigmax(image,4);
[scale1m, scale2m, scale3m, LL] = a2m(noisedimage);

descale1m = zeros(256,256,3);
descale2m = zeros(128,128,3);
descale3m = zeros(64,64,3);

switch method
    case {1}
        for i = 1:3
            for j = 1:256
                for k = 1:256
                    tempT = Tbo(sigma,sigmax1m(j,k,i));
                    descale1m(j,k,i) = wthresh(scale1m(j,k,i),'s',tempT);
                end
            end
        end
        for i = 1:3
            for j = 1:128
                for k = 1:128
                    tempT = Tbo(sigma,sigmax2m(j,k,i));
                    descale2m(j,k,i) = wthresh(scale2m(j,k,i),'s',tempT);
                end
            end
        end
        for i = 1:3
            for j = 1:64
                for k = 1:64
                    tempT = Tbo(sigma,sigmax3m(j,k,i));
                    descale3m(j,k,i) = wthresh(scale3m(j,k,i),'s',tempT);
                end
            end
        end
        
    case {2}
        for i = 1:3
            for j = 1:256
                for k = 1:256
                    tempT = Tb10(sigma,sigmax1m(j,k,i));
                    descale1m(j,k,i) = wthresh(scale1m(j,k,i),'s',tempT);
                end
            end
        end
        for i = 1:3
            for j = 1:128
                for k = 1:128
                    tempT = Tb10(sigma,sigmax2m(j,k,i));
                    descale2m(j,k,i) = wthresh(scale2m(j,k,i),'s',tempT);
                end
            end
        end
        for i = 1:3
            for j = 1:64
                for k = 1:64
                    tempT = Tb05(sigma,sigmax3m(j,k,i));
                    descale3m(j,k,i) = wthresh(scale3m(j,k,i),'s',tempT);
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
function [T] = Tbo(sigma,sigmax)

T = sigma^2/sigmax;

end

% Gaussian
function [T] = Tb20(sigma,sigmax)

x = sigmax/sigma;
T = sigma^2/sigmax*(0.0944/(x+0.0119)+0.78);

end

% Laplacian
function [T] = Tb10(sigma,sigmax)

x = sigmax/sigma;
T = sigma^2/sigmax*(-0.01903*x.^2+0.1642*x+0.8104);

end

% Beta = 0.5
function [T] = Tb05(sigma,sigmax)

x = sigmax/sigma;
T = sigma^2/sigmax*(0.3228*(x)+0.7691);

end