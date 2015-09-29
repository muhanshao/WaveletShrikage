function [ denoisedImage ] = band_Denoise( image, noisedimage, sigma, method, thrtype )
%% band_Denoise Denoise a image using wavelet thresholding method and calculate sigmax locally.
%   [denoisedimage] = band_Denoise(image, noisedimage, sigma, method, thrtype);
%   image: clean image
%   noisedimage: noised image
%   sigma: noise standard deviation
%   method = 1: calculate real sigmax using clean image, thresholding by original formula
%   method = 2: calculate real sigmax using clean image, thresholding by modified formula
%   method = 3: estimate sigmax using noised image, thrsholding by original formula
%   method = 4: estimate sigmax using noised image, thrsholding by modified formula
%   thrtype: thresholding type. 'h': hard-thresholding. 's':
%   soft-thresholding.

S = [64 64;64 64; 128 128;256 256;512 512];
[scale1, scale2, scale3] = dwt2_3level(noisedimage);

switch method
    case {1}
        [scale1r, scale2r, scale3r] = dwt2_3level(image);
        % Calculate sigmax
        sigmaxrS1HL = sqrt(sum(scale1r(1,:).^2)/65536);
        sigmaxrS1LH = sqrt(sum(scale1r(2,:).^2)/65536);
        sigmaxrS1HH = sqrt(sum(scale1r(3,:).^2)/65536);
        sigmaxrS2HL = sqrt(sum(scale2r(1,:).^2)/16384);
        sigmaxrS2LH = sqrt(sum(scale2r(2,:).^2)/16384);
        sigmaxrS2HH = sqrt(sum(scale2r(3,:).^2)/16384);
        sigmaxrS3HL = sqrt(sum(scale3r(2,:).^2)/4096);
        sigmaxrS3LH = sqrt(sum(scale3r(3,:).^2)/4096);
        sigmaxrS3HH = sqrt(sum(scale3r(4,:).^2)/4096);
        % Calculate Tb
        TbS1HL = Tbo(sigma,sigmaxrS1HL);
        TbS1LH = Tbo(sigma,sigmaxrS1LH);
        TbS1HH = Tbo(sigma,sigmaxrS1HH);
        TbS2HL = Tbo(sigma,sigmaxrS2HL);
        TbS2LH = Tbo(sigma,sigmaxrS2LH);
        TbS2HH = Tbo(sigma,sigmaxrS2HH);
        TbS3HL = Tbo(sigma,sigmaxrS3HL);
        TbS3LH = Tbo(sigma,sigmaxrS3LH);
        TbS3HH = Tbo(sigma,sigmaxrS3HH);
        % Thresholding
         % scale1
        scale11Thr = wthresh(scale1(1,:),thrtype,TbS1HL);
        scale12Thr = wthresh(scale1(2,:),thrtype,TbS1LH);
        scale13Thr = wthresh(scale1(3,:),thrtype,TbS1HH);
        % scale2
        scale21Thr = wthresh(scale2(1,:),thrtype,TbS2HL);
        scale22Thr = wthresh(scale2(2,:),thrtype,TbS2LH);
        scale23Thr = wthresh(scale2(3,:),thrtype,TbS2HH);
        % scale3
        scale32Thr = wthresh(scale3(2,:),thrtype,TbS3HL);
        scale33Thr = wthresh(scale3(3,:),thrtype,TbS3LH);
        scale34Thr = wthresh(scale3(4,:),thrtype,TbS3HH); 
        
    case {2}
        [scale1r, scale2r, scale3r] = dwt2_3level(image);
        % Calculate sigmax
        sigmaxrS1HL = sqrt(sum(scale1r(1,:).^2)/65536);
        sigmaxrS1LH = sqrt(sum(scale1r(2,:).^2)/65536);
        sigmaxrS1HH = sqrt(sum(scale1r(3,:).^2)/65536);
        sigmaxrS2HL = sqrt(sum(scale2r(1,:).^2)/16384);
        sigmaxrS2LH = sqrt(sum(scale2r(2,:).^2)/16384);
        sigmaxrS2HH = sqrt(sum(scale2r(3,:).^2)/16384);
        sigmaxrS3HL = sqrt(sum(scale3r(2,:).^2)/4096);
        sigmaxrS3LH = sqrt(sum(scale3r(3,:).^2)/4096);
        sigmaxrS3HH = sqrt(sum(scale3r(4,:).^2)/4096);
        % Calculate Tb
        TbS1HL = Tb051(sigma,sigmaxrS1HL);
        TbS1LH = Tb051(sigma,sigmaxrS1LH);
        TbS1HH = Tb051(sigma,sigmaxrS1HH);
        TbS2HL = Tb051(sigma,sigmaxrS2HL);
        TbS2LH = Tb051(sigma,sigmaxrS2LH);
        TbS2HH = Tb051(sigma,sigmaxrS2HH);
        TbS3HL = Tb051(sigma,sigmaxrS3HL);
        TbS3LH = Tb051(sigma,sigmaxrS3LH);
        TbS3HH = Tb051(sigma,sigmaxrS3HH);
        % Thresholding
         % scale1
        scale11Thr = wthresh(scale1(1,:),thrtype,TbS1HL);
        scale12Thr = wthresh(scale1(2,:),thrtype,TbS1LH);
        scale13Thr = wthresh(scale1(3,:),thrtype,TbS1HH);
        % scale2
        scale21Thr = wthresh(scale2(1,:),thrtype,TbS2HL);
        scale22Thr = wthresh(scale2(2,:),thrtype,TbS2LH);
        scale23Thr = wthresh(scale2(3,:),thrtype,TbS2HH);
        % scale3
        scale32Thr = wthresh(scale3(2,:),thrtype,TbS3HL);
        scale33Thr = wthresh(scale3(3,:),thrtype,TbS3LH);
        scale34Thr = wthresh(scale3(4,:),thrtype,TbS3HH); 
    
    case {3}
        % Calculate sigmayhat^2
        % scale1
        sigmayhatS1HL = sum(scale1(1,:).^2)/65536;
        sigmayhatS1LH = sum(scale1(2,:).^2)/65536;
        sigmayhatS1HH = sum(scale1(3,:).^2)/65536;
        % scale2
        sigmayhatS2HL = sum(scale2(1,:).^2)/16384;
        sigmayhatS2LH = sum(scale2(2,:).^2)/16384;
        sigmayhatS2HH = sum(scale2(3,:).^2)/16384;
        % scale3
        sigmayhatS3HL = sum(scale3(2,:).^2)/4096;
        sigmayhatS3LH = sum(scale3(3,:).^2)/4096;
        sigmayhatS3HH = sum(scale3(4,:).^2)/4096;

        % Calculate sigmaxhat
        % scale1
        sigmaxhatS1HL = sqrt(max(sigmayhatS1HL-sigma^2,0));
        sigmaxhatS1LH = sqrt(max(sigmayhatS1LH-sigma^2,0));
        sigmaxhatS1HH = sqrt(max(sigmayhatS1HH-sigma^2,0));
        % scale2
        sigmaxhatS2HL = sqrt(max(sigmayhatS2HL-sigma^2,0));
        sigmaxhatS2LH = sqrt(max(sigmayhatS2LH-sigma^2,0));
        sigmaxhatS2HH = sqrt(max(sigmayhatS2HH-sigma^2,0));
        % scale3
        sigmaxhatS3HL = sqrt(max(sigmayhatS3HL-sigma^2,0));
        sigmaxhatS3LH = sqrt(max(sigmayhatS3LH-sigma^2,0));
        sigmaxhatS3HH = sqrt(max(sigmayhatS3HH-sigma^2,0));

        % Calculate Tbhat
        % scale1
        if sigmaxhatS1HL == 0
            TbhatS1HL = max(scale1(1,:));
        else
            TbhatS1HL = Tbo(sigma,sigmaxhatS1HL);
        end
        if sigmaxhatS1LH == 0
            TbhatS1LH = max(scale1(2,:));
        else
            TbhatS1LH = Tbo(sigma,sigmaxhatS1LH);
        end
        if sigmaxhatS1HH == 0
            TbhatS1HH = max(scale1(3,:));
        else
            TbhatS1HH = Tbo(sigma,sigmaxhatS1HH);
        end
        % scale2
        if sigmaxhatS2HL == 0
            TbhatS2HL = max(scale2(1,:));
        else
            TbhatS2HL = Tbo(sigma,sigmaxhatS2HL);
        end
        if sigmaxhatS2LH == 0
            TbhatS2LH = max(scale2(2,:));
        else
            TbhatS2LH = Tbo(sigma,sigmaxhatS2LH);
        end
        if sigmaxhatS2HH == 0
            TbhatS2HH = max(scale2(3,:));
        else
            TbhatS2HH = Tbo(sigma,sigmaxhatS2HH);
        end
        % scale3
        if sigmaxhatS3HL == 0
            TbhatS3HL = max(scale3(2,:));
        else
            TbhatS3HL = Tbo(sigma,sigmaxhatS3HL);
        end
        if sigmaxhatS3LH == 0
            TbhatS3LH = max(scale3(3,:));
        else
            TbhatS3LH = Tbo(sigma,sigmaxhatS3LH);
        end
        if sigmaxhatS3HH == 0
            TbhatS3HH = max(scale3(4,:));
        else
            TbhatS3HH = Tbo(sigma,sigmaxhatS3HH);
        end

        % Thresholding
        % scale1
        scale11Thr = wthresh(scale1(1,:),thrtype,TbhatS1HL);
        scale12Thr = wthresh(scale1(2,:),thrtype,TbhatS1LH);
        scale13Thr = wthresh(scale1(3,:),thrtype,TbhatS1HH);
        % scale2
        scale21Thr = wthresh(scale2(1,:),thrtype,TbhatS2HL);
        scale22Thr = wthresh(scale2(2,:),thrtype,TbhatS2LH);
        scale23Thr = wthresh(scale2(3,:),thrtype,TbhatS2HH);
        % scale3
        scale32Thr = wthresh(scale3(2,:),thrtype,TbhatS3HL);
        scale33Thr = wthresh(scale3(3,:),thrtype,TbhatS3LH);
        scale34Thr = wthresh(scale3(4,:),thrtype,TbhatS3HH); 
        
    case {4}
        % Calculate sigmayhat^2
        % scale1
        sigmayhatS1HL = sum(scale1(1,:).^2)/65536;
        sigmayhatS1LH = sum(scale1(2,:).^2)/65536;
        sigmayhatS1HH = sum(scale1(3,:).^2)/65536;
        % scale2
        sigmayhatS2HL = sum(scale2(1,:).^2)/16384;
        sigmayhatS2LH = sum(scale2(2,:).^2)/16384;
        sigmayhatS2HH = sum(scale2(3,:).^2)/16384;
        % scale3
        sigmayhatS3HL = sum(scale3(2,:).^2)/4096;
        sigmayhatS3LH = sum(scale3(3,:).^2)/4096;
        sigmayhatS3HH = sum(scale3(4,:).^2)/4096;

        % Calculate sigmaxhat
        % scale1
        sigmaxhatS1HL = sqrt(max(sigmayhatS1HL-sigma^2,0));
        sigmaxhatS1LH = sqrt(max(sigmayhatS1LH-sigma^2,0));
        sigmaxhatS1HH = sqrt(max(sigmayhatS1HH-sigma^2,0));
        % scale2
        sigmaxhatS2HL = sqrt(max(sigmayhatS2HL-sigma^2,0));
        sigmaxhatS2LH = sqrt(max(sigmayhatS2LH-sigma^2,0));
        sigmaxhatS2HH = sqrt(max(sigmayhatS2HH-sigma^2,0));
        % scale3
        sigmaxhatS3HL = sqrt(max(sigmayhatS3HL-sigma^2,0));
        sigmaxhatS3LH = sqrt(max(sigmayhatS3LH-sigma^2,0));
        sigmaxhatS3HH = sqrt(max(sigmayhatS3HH-sigma^2,0));

        % Calculate Tbhat
        % scale1
        if sigmaxhatS1HL == 0
            TbhatS1HL = max(scale1(1,:));
        else
            TbhatS1HL = Tb051(sigma,sigmaxhatS1HL);
        end
        if sigmaxhatS1LH == 0
            TbhatS1LH = max(scale1(2,:));
        else
            TbhatS1LH = Tb051(sigma,sigmaxhatS1LH);
        end
        if sigmaxhatS1HH == 0
            TbhatS1HH = max(scale1(3,:));
        else
            TbhatS1HH = Tb051(sigma,sigmaxhatS1HH);
        end
        % scale2
        if sigmaxhatS2HL == 0
            TbhatS2HL = max(scale2(1,:));
        else
            TbhatS2HL = Tb051(sigma,sigmaxhatS2HL);
        end
        if sigmaxhatS2LH == 0
            TbhatS2LH = max(scale2(2,:));
        else
            TbhatS2LH = Tb051(sigma,sigmaxhatS2LH);
        end
        if sigmaxhatS2HH == 0
            TbhatS2HH = max(scale2(3,:));
        else
            TbhatS2HH = Tb051(sigma,sigmaxhatS2HH);
        end
        % scale3
        if sigmaxhatS3HL == 0
            TbhatS3HL = max(scale3(2,:));
        else
            TbhatS3HL = Tb051(sigma,sigmaxhatS3HL);
        end
        if sigmaxhatS3LH == 0
            TbhatS3LH = max(scale3(3,:));
        else
            TbhatS3LH = Tb051(sigma,sigmaxhatS3LH);
        end
        if sigmaxhatS3HH == 0
            TbhatS3HH = max(scale3(4,:));
        else
            TbhatS3HH = Tb051(sigma,sigmaxhatS3HH);
        end

        % Thresholding
        % scale1
        scale11Thr = wthresh(scale1(1,:),thrtype,TbhatS1HL);
        scale12Thr = wthresh(scale1(2,:),thrtype,TbhatS1LH);
        scale13Thr = wthresh(scale1(3,:),thrtype,TbhatS1HH);
        % scale2
        scale21Thr = wthresh(scale2(1,:),thrtype,TbhatS2HL);
        scale22Thr = wthresh(scale2(2,:),thrtype,TbhatS2LH);
        scale23Thr = wthresh(scale2(3,:),thrtype,TbhatS2HH);
        % scale3
        scale32Thr = wthresh(scale3(2,:),thrtype,TbhatS3HL);
        scale33Thr = wthresh(scale3(3,:),thrtype,TbhatS3LH);
        scale34Thr = wthresh(scale3(4,:),thrtype,TbhatS3HH); 
            
    otherwise
        disp('Unknown method!');     
end

CTh = [scale3(1,:) scale32Thr scale33Thr scale34Thr scale21Thr scale22Thr scale23Thr scale11Thr scale12Thr scale13Thr];
denoisedImage = waverec2(CTh,S,'haar');

end


function [T] = Tbo(sigma,sigmax)

T = sigma^2/sigmax;

end

function [T] = Tb051(sigma,sigmax)

x = sigmax/sigma;
T = sigma^2/sigmax*(-0.06958*x^2+0.6081*x+0.5672);

end

