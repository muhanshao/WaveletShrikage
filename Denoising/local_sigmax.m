function [ sigmax1m, sigmax2m, sigmax3m, flag1m, flag2m, flag3m ] = local_sigmax( image,windowsize,sigma,type )
%% local_sigmax Calculating local sigmax
%   [ sigmax1m, sigmax2m, sigmax3m ] = local_cal( image,windowsize,type )
%   sigmax1~3m are 3-D matrix containing every local sigmax.
%   windowsize typically is 3 or 4.
%   sigma: noise level.
%   type = 1, clean image. type = 2, noised image.

sigmay1m = zeros(256,256,3);
sigmay2m = zeros(128,128,3);
sigmay3m = zeros(64,64,3);

sigmax1m = zeros(256,256,3);
sigmax2m = zeros(128,128,3);
sigmax3m = zeros(64,64,3);

flag1m = zeros(256,256,3);
flag2m = zeros(128,128,3);
flag3m = zeros(64,64,3);

[ scale1m, scale2m, scale3m, LL ] = a2m( image );

switch type
    case {1}
        % Calculate level 1
        for i = 1:3
            for j = 1:256
                for k = 1:256
                   if j-windowsize<1 
                       if k-windowsize<1
                           tempx = scale1m(1:j+windowsize,1:k+windowsize,i);
                       elseif k+windowsize>256
                           tempx = scale1m(1:j+windowsize,k-windowsize:256,i);
                       else
                           tempx = scale1m(1:j+windowsize,k-windowsize:k+windowsize,i);
                       end
                   elseif j+windowsize>256
                       if k-windowsize<1
                           tempx = scale1m(j-windowsize:256,1:k+windowsize,i);
                       elseif k+windowsize>256
                           tempx = scale1m(j-windowsize:256,k-windowsize:256,i);
                       else
                           tempx = scale1m(j-windowsize:256,k-windowsize:k+windowsize,i);
                       end
                   else
                       if k-windowsize<1
                           tempx = scale1m(j-windowsize:j+windowsize,1:k+windowsize,i);
                       elseif k+windowsize>256
                           tempx = scale1m(j-windowsize:j+windowsize,k-windowsize:256,i);
                       else
                           tempx = scale1m(j-windowsize:j+windowsize,k-windowsize:k+windowsize,i);
                       end
                   end
                   sigmax1m(j,k,i) = Calsigmax(tempx);
                end
            end
        end

        % Calculate level 2
        for i = 1:3
            for j = 1:128
                for k = 1:128
                   if j-windowsize<1 
                       if k-windowsize<1
                           tempx = scale2m(1:j+windowsize,1:k+windowsize,i);
                       elseif k+windowsize>128
                           tempx = scale2m(1:j+windowsize,k-windowsize:128,i);
                       else
                           tempx = scale2m(1:j+windowsize,k-windowsize:k+windowsize,i);
                       end
                   elseif j+windowsize>128
                       if k-windowsize<1
                           tempx = scale2m(j-windowsize:128,1:k+windowsize,i);
                       elseif k+windowsize>128
                           tempx = scale2m(j-windowsize:128,k-windowsize:128,i);
                       else
                           tempx = scale2m(j-windowsize:128,k-windowsize:k+windowsize,i);
                       end
                   else
                       if k-windowsize<1
                           tempx = scale2m(j-windowsize:j+windowsize,1:k+windowsize,i);
                       elseif k+windowsize>128
                           tempx = scale2m(j-windowsize:j+windowsize,k-windowsize:128,i);
                       else
                           tempx = scale2m(j-windowsize:j+windowsize,k-windowsize:k+windowsize,i);
                       end
                   end
                   sigmax2m(j,k,i) = Calsigmax(tempx);
                end
            end
        end

        % Calculate level 3
        for i = 1:3
            for j = 1:64
                for k = 1:64
                   if j-windowsize<1 
                       if k-windowsize<1
                           tempx = scale3m(1:j+windowsize,1:k+windowsize,i);
                       elseif k+windowsize>64
                           tempx = scale3m(1:j+windowsize,k-windowsize:64,i);
                       else
                           tempx = scale3m(1:j+windowsize,k-windowsize:k+windowsize,i);
                       end
                   elseif j+windowsize>64
                       if k-windowsize<1
                           tempx = scale3m(j-windowsize:64,1:k+windowsize,i);
                       elseif k+windowsize>64
                           tempx = scale3m(j-windowsize:64,k-windowsize:64,i);
                       else
                           tempx = scale3m(j-windowsize:64,k-windowsize:k+windowsize,i);
                       end
                   else
                       if k-windowsize<1
                           tempx = scale3m(j-windowsize:j+windowsize,1:k+windowsize,i);
                       elseif k+windowsize>64
                           tempx = scale3m(j-windowsize:j+windowsize,k-windowsize:64,i);
                       else
                           tempx = scale3m(j-windowsize:j+windowsize,k-windowsize:k+windowsize,i);
                       end
                   end
                   sigmax3m(j,k,i) = Calsigmax(tempx);
                end
            end
        end

    case {2}
        % Calculate level 1
        for i = 1:3
            for j = 1:256
                for k = 1:256
                   if j-windowsize<1 
                       if k-windowsize<1
                           tempy = scale1m(1:j+windowsize,1:k+windowsize,i);
                       elseif k+windowsize>256
                           tempy = scale1m(1:j+windowsize,k-windowsize:256,i);
                       else
                           tempy = scale1m(1:j+windowsize,k-windowsize:k+windowsize,i);
                       end
                   elseif j+windowsize>256
                       if k-windowsize<1
                           tempy = scale1m(j-windowsize:256,1:k+windowsize,i);
                       elseif k+windowsize>256
                           tempy = scale1m(j-windowsize:256,k-windowsize:256,i);
                       else
                           tempy = scale1m(j-windowsize:256,k-windowsize:k+windowsize,i);
                       end
                   else
                       if k-windowsize<1
                           tempy = scale1m(j-windowsize:j+windowsize,1:k+windowsize,i);
                       elseif k+windowsize>256
                           tempy = scale1m(j-windowsize:j+windowsize,k-windowsize:256,i);
                       else
                           tempy = scale1m(j-windowsize:j+windowsize,k-windowsize:k+windowsize,i);
                       end
                   end
                   sigmay1m(j,k,i) = Calsigmax(tempy);
                   sigmax1m(j,k,i) = sqrt(max(sigmay1m(j,k,i)^2-sigma^2,0));
                   if sigmax1m(j,k,i) == 0
                      flag1m(j,k,i) = 1;
                      sigmax1m(j,k,i) = max(tempy(:));
                   end
                end
            end
        end

        % Calculate level 2
        for i = 1:3
            for j = 1:128
                for k = 1:128
                   if j-windowsize<1 
                       if k-windowsize<1
                           tempy = scale2m(1:j+windowsize,1:k+windowsize,i);
                       elseif k+windowsize>128
                           tempy = scale2m(1:j+windowsize,k-windowsize:128,i);
                       else
                           tempy = scale2m(1:j+windowsize,k-windowsize:k+windowsize,i);
                       end
                   elseif j+windowsize>128
                       if k-windowsize<1
                           tempy = scale2m(j-windowsize:128,1:k+windowsize,i);
                       elseif k+windowsize>128
                           tempy = scale2m(j-windowsize:128,k-windowsize:128,i);
                       else
                           tempy = scale2m(j-windowsize:128,k-windowsize:k+windowsize,i);
                       end
                   else
                       if k-windowsize<1
                           tempy = scale2m(j-windowsize:j+windowsize,1:k+windowsize,i);
                       elseif k+windowsize>128
                           tempy = scale2m(j-windowsize:j+windowsize,k-windowsize:128,i);
                       else
                           tempy = scale2m(j-windowsize:j+windowsize,k-windowsize:k+windowsize,i);
                       end
                   end
                   sigmay2m(j,k,i) = Calsigmax(tempy);
                   sigmax2m(j,k,i) = sqrt(max(sigmay2m(j,k,i)^2-sigma^2,0));
                   if sigmax2m(j,k,i) == 0
                      flag2m(j,k,i) = 1;
                      sigmax2m(j,k,i) = max(tempy(:));
                   end
                end
            end
        end

        % Calculate level 3
        for i = 1:3
            for j = 1:64
                for k = 1:64
                   if j-windowsize<1 
                       if k-windowsize<1
                           tempy = scale3m(1:j+windowsize,1:k+windowsize,i);
                       elseif k+windowsize>64
                           tempy = scale3m(1:j+windowsize,k-windowsize:64,i);
                       else
                           tempy = scale3m(1:j+windowsize,k-windowsize:k+windowsize,i);
                       end
                   elseif j+windowsize>64
                       if k-windowsize<1
                           tempy = scale3m(j-windowsize:64,1:k+windowsize,i);
                       elseif k+windowsize>64
                           tempy = scale3m(j-windowsize:64,k-windowsize:64,i);
                       else
                           tempy = scale3m(j-windowsize:64,k-windowsize:k+windowsize,i);
                       end
                   else
                       if k-windowsize<1
                           tempy = scale3m(j-windowsize:j+windowsize,1:k+windowsize,i);
                       elseif k+windowsize>64
                           tempy = scale3m(j-windowsize:j+windowsize,k-windowsize:64,i);
                       else
                           tempy = scale3m(j-windowsize:j+windowsize,k-windowsize:k+windowsize,i);
                       end
                   end
                   sigmay3m(j,k,i) = Calsigmax(tempy);
                   sigmax3m(j,k,i) = sqrt(max(sigmay3m(j,k,i)^2-sigma^2,0));
                   if sigmax3m(j,k,i) == 0
                      flag3m(j,k,i) = 1;
                      sigmax3m(j,k,i) = max(tempy(:));
                   end
                end
            end
        end
        
    otherwise
        disp('Unknown method!');        
end

end

function [ sigmax ] = Calsigmax( x )

x = x(:);
N = length(x);
sigmax = sqrt(sum(x.^2)/N);

end


