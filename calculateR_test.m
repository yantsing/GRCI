% Give the channels H, noise variance \sigma^2, and regularization
% parameter \alpha, calculate the true secrecy sum-rate using RCI precoding. 


function [R,R1] = calculateR_test(H, sigma2, alpha)
UTEST = 1;
DEBUG_A_B = 0
K = size(H,1);
N = size(H,2);
I_K = eye(K);
I_N = eye(N);
W = H' * inv(H*H' + alpha * I_K);
zeta = real(trace(W' * W))
gamm = real(trace(H'*H*(H'*H + alpha * I_N)^(-2)))

X = (H'*H + alpha * I_N)^(-1);

for k = 1 : K
    h_k = H(k,:);
    
    if k == 1
        H_k = H(2:K, :);
    elseif k == K
        H_k = H(1:K-1, :);
    else
        index = [1:k-1,k+1:K]
        H_k = H(index,:);
    end
    
    temp = real(h_k* X * H_k' * H_k * X* h_k');
    
    SINR = (abs(h_k * X * h_k'))^2/(zeta * sigma2 + temp);
    R(k) = log2(1 + SINR);
end



Xi = alpha * eye(K);

[U,D,V] = svd(H);
Lambda = D*D;
e = eye(K,K);
for k = 1 : k
    I(:,:,k) = e(:,k)*e(:,k)';
end

W = H'*inv(H*H'+ U*Xi*U');
%%%%%%%%%%Calculating vector x %%%%%%%%%%%%%%%%%%%%%%
for k = 1 : K
    x(k) = 1/(Lambda(k,k)+Xi(k,k));
end
x = x';

%%%%%%%%%%Calculating vector b %%%%%%%%%%%%%%%%%%%%%%
b = diag(D);

%%%%%%%%%%Calculating B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = b*b';

%%%%%%%%%%Calculating vectors a_kj %%%%%%%%%%%%%%%%%%%%%
for k = 1 : K
    for j = 1 : K
        a(:,k,j) = (H(k,:) * V * D * diag(U(j,:)'))';
    end
end

%%%%%%%%%%%%%%%Claculating t_k%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : K
    if UTEST == 1
        % Teh following quantities should be the same.
        t_test(k)= H(k,:)*H'*inv(H*H'+ U*Xi*U')*e(:,k)
        t_test1(k) = H(k,:)* W(:,k)
        t(k) = H(k,:) * H'* U * inv(Lambda + Xi)*U(k,:)'
        t(k) = H(k,:) * H'* U * diag(U(k,:)') * x
        
        t(k) = H(k,:) * V * D * diag(U(k,:)')* x
        t(k) = a(:,k,k)'*x
    else
        t(k) = a(:,k,k)'*x;
    end
end

%%%%%%%%%%%%%%%%Claculating t_{\tilde{k}} %%%%%%%%%%%%%%%%%%%
for k = 1 : K
    temp_test = 0;
    temp = 0;
    temp_test1 = 0;
    temp_test2 = 0;
    for j = 1 : K
        if j ~= k
            if UTEST == 1
                temp_test = temp_test + abs(H(k,:)* W(:,j))^2;
                temp_test1 = temp_test1 + abs(H(k,:) * V * D * diag(U(j,:)')* x)^2;
                temp = temp + abs(a(:,k,j)'*x)^2;
            else
                temp = temp + abs(a(:,k,j)'*x)^2;
            end
        end
    end
    tTilde2_test(k) = temp_test
    tTilde2_test(k) = temp_test1
    tTilde2(k)=temp
end

%%%%%%%%%%%%%%%%%Calculating A_k,A_{\tilde{k}}%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : K
    A(:,:,k) = a(:,k,k)* a(:,k,k)';
    temp = zeros(K);
    for j = 1 : K
        if j ~= k
                temp = temp + a(:,k,j)*a(:,k,j)';
        end
    end
    ATilde(:,:,k) = temp
end

%%%%%%%%%%%%%%%Calculating SINR_k%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    R1(k) = log2(1 + numerator/denomenator);
end
end




