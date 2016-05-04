% Give the channels H, noise variance \sigma^2, and regularization
% parameter \alpha, calculate the true secrecy sum-rate using RCI precoding. 


function R = sumCapacity(H, sigma2)
% P = 1;
% [U, Lam, V] = svd(H'* H);
% 
% LamInv = 1./diag(Lam);
% wline = wfill(LamInv, P, 0.001);
% lam_Q = wline - LamInv;
% lam_Q(lam_Q < 0) = 0;
% Q = U*diag(lam_Q)*V';

TransmitPower = 1;
N = size(H,2);

K = size(H,1);

cvx_begin quiet
    
variable covQ(N,N) hermitian
maximize(log_det(eye(K) + sigma2^(-1)*H*covQ*H'));
subject to
covQ == hermitian_semidefinite(N);
abs(trace(covQ)) <= TransmitPower;

cvx_end

Q = covQ;

R = log2(det(eye(K) + sigma2^(-1)*H*Q*H'));


