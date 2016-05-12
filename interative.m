function [ xi ] = interative( lambda,alpha)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

K = length(lambda);

y = ones(1,K);
for l = 1 : K
    y(l) = lambda(l)/(lambda(l) + alpha);
end
sigma2 = 0.1;
sigma2K = sigma2 * K;
for it = 1 :100
for k = 1 :K
    tempN = 0;
    tempD = 0;
    for i = 1 : K
        if i ~= k
            tempN = sigma2K * y(i)^2/lambda(i) + y(i)^2  + tempN;
            tempD = sigma2K * y(i)/lambda(k) + y(i) + tempD;
        end
    end
    tempy(k) = tempN/tempD;
end
y = tempy;
end

y = y/max(y);

for k = 1 : K
    xi(k) =  (1/y(k) - 1) * lambda(k);
end
end

