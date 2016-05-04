%H is channele, W is the weight vector%

function R = calculateR(x, sigma2, A, ATilde,tTilde2,Lambda,a, UTEST)
K = size(A, 1);
for k = 1 : K
    if UTEST == 1 
        same11 =  norm(a(:,k,k)'*x)^2
        same12  = x'*A(:,:,k)*x
        same13 = (real(a(:,k,k))'*x)^2
        
        same21 = tTilde2(k) + sigma2 * trace(Lambda * x * x')
        same22 =  x' * (ATilde(:,:,k) + sigma2 * Lambda) * x
        same23 = norm(sqrtm(ATilde(:,:,k) + sigma2 * Lambda) * x)^2
        numerator = same13;
        denomenator = same23;
        
    else
        numerator = real(x'*A(:,:,k)*x);
        denomenator = real(x' * (ATilde(:,:,k) + sigma2 * Lambda) * x);
    end
    R(k) = log2(1 + numerator/denomenator);
end
end