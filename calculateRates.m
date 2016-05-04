% Give the channels H, weight vectors W

function R = calculateRates(H,W,sigma2)
K = size(H,1);
for k = 1 : K
    desired = abs(H(k,:)* W(:,k))^2;
    interference = 0;
    for j = 1 : K
        if j ~= k
            interference = interference + abs(H(k,:)*W(:,j))^2;
        end
    end
    R(k) = log2(1 + desired / (interference + sigma2));
end

end




