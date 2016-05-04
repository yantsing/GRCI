H = rand(4,4);
I = eye(4);
DEBUT_A_B = 1
alpha = 0.01;
sigma = 0.1;
W = H' * inv(H*H' + alpha * I);
zeta = real(trace(W' * W))

for k = 1 : K
   if k == 1
        H_k(:,:,k) = H(2:K, :);
    elseif k == K
        H_k (:,:,k)= H(1:K-1, :);
    else
        index = [1:k-1,k+1:K];
        H_k(:,:,k) = H(index,:);
    end  
end

for k = 1 : K
    Z = H_k(:,:,k)' * H_k(:,:,k);
    I = size(Z, 1);
    
    
    
    A(k)= H(k,:) * inv(Z + alpha  * I) * H(k,:)';
    if DEBUT_A_B == 1
    test1 = H(k,:) * inv(H' * H + alpha * I) * H(k,:)'
    test2 = A(k)/ (1 + A(k))
    end
    
    test3 = real(H(k,:) * inv(H' * H + alpha * I) * H_k(:,:,k)' * H_k(:,:,k) * inv(H' * H + alpha * I)* H(k,:)')
    
    B(k) = H(k,:) * inv(Z + alpha  * I) * Z* inv(Z + alpha  * I) * H(k,:)';
    if DEBUT_A_B == 1
        test3 = real(H(k,:) * inv(H' * H + alpha * I) * H_k(:,:,k)' * H_k(:,:,k) * inv(H' * H + alpha * I)* H(k,:)')
        test4 =  B(k)/(1+ A(k))^2
    end
    SINR1(k) = A(k)^2/(B(k) + zeta * sigma2 * (1 + A(k))^2);
    SINR1_tilde(k) = B(k)/(zeta * sigma2 * (1 + A(k))^2);
end

R1 = 0;
for i = 1 : k
    R1 = R1 + max(log2(1 + SINR1(k))-log2(1 + SINR1_tilde(k)),0)
end
