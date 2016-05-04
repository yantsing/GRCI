% Give the channels H, weight vectors W

function R = calculateCapacity(H,W,sigma2)
K = size(H,1);
R = log2(det(eye(K) + sigma2^(-1)*H*W*W'*H'));
end




