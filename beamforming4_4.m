% Give the channels H, noise variance \sigma^2, and regularization
% parameter \alpha, calculate the rates. 
function result = beamforming4_4(H, sigma2)
K = size(H,1);
N = size(H,2);

t_max = 10;
t_min = 0;
threshold = 0.001
result = zeros(N,K);
for k = 1 : K
    Hk(:,:,k) = H(k,:)'*H(k,:);
end



while t_max - t_min > threshold

t = (t_max + t_min)/2;
cvx_begin
cvx_precision best
variable  W1(N,N) complex semidefinite
variable  W2(N,N) complex semidefinite
variable  W3(N,N) complex semidefinite
variable  W4(N,N) complex semidefinite



minimize(trace(W1) + trace(W2) + trace(W3) + trace(W4))

subject to

trace(Hk(:,:,1) * W1)  -  t^2 * (trace(Hk(:,:,1) * W2) + trace(Hk(:,:,1) * W3)+trace(Hk(:,:,1) * W4)) >= sigma2 * t^2;
trace(Hk(:,:,2) * W2)  -  t^2 * (trace(Hk(:,:,2) * W1) + trace(Hk(:,:,2) * W3)+trace(Hk(:,:,2) * W4)) >= sigma2 * t^2;
trace(Hk(:,:,3) * W3)  -  t^2 * (trace(Hk(:,:,3) * W1) + trace(Hk(:,:,3) * W2)+trace(Hk(:,:,3) * W4)) >= sigma2 * t^2;
trace(Hk(:,:,4) * W4)  -  t^2 * (trace(Hk(:,:,4) * W1) + trace(Hk(:,:,4) * W2)+trace(Hk(:,:,4) * W3)) >= sigma2 * t^2;

trace(W1) + trace(W2) + trace(W3) + trace(W4) <= 1

cvx_end
if strcmp(cvx_status, 'Infeasible') == 1 || strcmp(cvx_status, 'Failed') 
    t_max = t;
else
    t_min = t;
    rW1 = W1;
    rW2 = W2;
    rW3 = W3;
    rW4 = W4;
end
end
[w1,d1]= svd(rW1);
[w2,d2]= svd(rW2);
[w3,d3]= svd(rW3);
[w4,d4]= svd(rW4);

result = [w1(:,1),w2(:,1),w3(:,1),w4(:,1)];
end