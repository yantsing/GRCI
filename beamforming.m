% Give the channels H, noise variance \sigma^2, and regularization
% parameter \alpha, calculate the rates. 
function W = beamforming4_4(H, sigma2)
K = size(H,1);
N = size(H,2);

W = zeros(N,K);
t_max = 10;
t_min = 0;
threshold = 0.001

for k = 1 : K
    Hk(:,:,k) = H(k,:)'*H(k,:);
end



while t_max - t_min > threshold

t = (t_max + t_min)/2;
cvx_begin quiet
cvx_precision best
variable  w1(K) complex
variable  w2(K) complex
variable  w3(K) complex
variable  w4(K) complex



minimize(w1'*w1 + w2' * w2 + w3' * w3 + w4' * w4)

subject to

real(H(1,:)*w1) >= t * norm([sigma2, H(1,:)*w2, H(1,:)*w2, H(1,:)*w3]) 

w1'*w1 + w2' * w2 + w3' * w3 + w4' * w4 <= 1

cvx_end
if strcmp(cvx_status, 'Infeasible') == 1 || strcmp(cvx_status, 'Failed') 
    t_max = t;
else
    t_min = t;
    W = [w1,w2,w3,w4];
end


end