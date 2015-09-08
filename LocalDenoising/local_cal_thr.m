function [ sigmax1m, sigmax2m, sigmax3m, beta1m, beta2m, beta3m ] = local_cal_thr( image,windowsize, thr, r,beta )
%local_cal Calculating local sigmax and beta
%   [ sigmax1m, sigmax2m, sigmax3m, beta1m, beta2m, beta3m ] = local_cal( image,windowsize, r,beta )
%   sigmax1~3m are 3-D matrix containing every local sigmax.
%   beta1~3m are 3-D matrix containing every local beta.
%   windowsize typically is 3 or 4.
%   r and beta are lookup table.

sigmax1m = zeros(256,256,3);
sigmax2m = zeros(128,128,3);
sigmax3m = zeros(64,64,3);
beta1m = zeros(256,256,3);
beta2m = zeros(128,128,3);
beta3m = zeros(64,64,3);

[ scale1m, scale2m, scale3m, LL ] = a2m( image );
scale1m = scale1m.*(scale1m>thr(1));
scale2m = scale2m.*(scale2m>thr(2));
scale3m = scale3m.*(scale3m>thr(3));

%% Calculate level 1
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
           beta1m(j,k,i) = Calbeta(tempx,r,beta);
        end
    end
end


%% Calculate level 2
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
           beta2m(j,k,i) = Calbeta(tempx,r,beta);
        end
    end
end

%% Calculate level 3
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
           beta3m(j,k,i) = Calbeta(tempx,r,beta);
        end
    end
end


end

function [ sigmax ] = Calsigmax( x )

x = x(:);
N = length(x);
sigmax = sqrt(sum(x.^2)/N);

end


