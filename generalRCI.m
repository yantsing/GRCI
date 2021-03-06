% This program simulate the CoMP conference paper.

clear all
DEBUG = 0
UTEST = 0
test_num = 10


t0 = clock;
%Parameter: # of MS, # of BS, # of receiving antenna, # of transmitting
%antenna.
N = 16; 
K = 4;
Ptr = 1;
beta = K/N;

SINR_SET = [-10]
sumC = zeros(length(SINR_SET),1);
sumR = zeros(length(SINR_SET),1);
sumRMax = zeros(length(SINR_SET),1);
sumRb = zeros(length(SINR_SET),1);
sumR_alpha = zeros(length(SINR_SET),1);
minR = zeros(length(SINR_SET),1);
minRMax = zeros(length(SINR_SET),1);
minRb = zeros(length(SINR_SET),1);
minR_alpha = zeros(length(SINR_SET),1);
sumAlpha = zeros(length(SINR_SET),1);
sumXi = zeros(length(SINR_SET),1);

for test = 1 :test_num
    test
for SINR = SINR_SET
sigma2 = 10^(-SINR/10);
rho = Ptr/sigma2;
alpha = beta/rho*N;


%%%%%%%%%%%%%%%%%%%Generate h(:,k)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1: K
%     l(k) = rand()
%     generatePathLoss
    l(k) = 1;
end
% l = generatePathLoss(K);
if DEBUG == 1
    l(1) = 1;
    l(2) = 0.1;
    l(3) = 0.1;
    l(4) = 0.1;
end

%Generate the channel h_k
for k = 1: K
    temp = (randn(1,N)+sqrt(-1)*randn(1,N))/sqrt(2);
    H(k,:) = l(k) * temp;
end

Xi = diag(rand(1,K));

%[U,D,V] = svd(H);
[U,Lambda] = eig(H*H');
e = eye(K,K);
for k = 1 : k
    I(:,:,k) = e(:,k)*e(:,k)';
end

W = H'*inv(H*H'+ U*Xi*U');
zeta = real(trace(W' * W));
W1  = W/sqrt(zeta);
%%%%%%%%%%Calculating vector x %%%%%%%%%%%%%%%%%%%%%%
for k = 1 : K
    x_temp(k) = 1/(Lambda(k,k)+Xi(k,k));
end
x = x_temp';

%%%%%%%%%%Calculating vector b %%%%%%%%%%%%%%%%%%%%%%
b = diag(Lambda);

%%%%%%%%%%Calculating B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = b*b';

%%%%%%%%%%Calculating vectors a_kj %%%%%%%%%%%%%%%%%%%%%
for k = 1 : K
    for j = 1 : K
        a(:,k,j) = (H(k,:) * H' * U * diag(U(j,:)'))';
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
                temp = temp + abs(a(:,k,j)'*x)^2;
            else
                temp = temp + abs(a(:,k,j)'*x)^2;
            end
        end
    end
    if UTEST == 1
        tTilde2_test(k) = temp_test
        tTilde2_test(k) = temp
    end
    tTilde2(k)=temp;
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
    ATilde(:,:,k) = temp;
end

%%%%%%%%%%%%%%%Calculating SINR_k%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : K
    if UTEST == 1 
        same11 =  norm(a(:,k,k)'*x)^2
        same12  = x'*A(:,:,k)*x
        
        same21 = tTilde2(k) + sigma2 * trace(Lambda * x * x')
        same22 =  x' * (ATilde(:,:,k) + sigma2 * Lambda) * x
        numerator = real(x'*A(:,:,k)*x)
        desired = abs(H(k,:)* W1(:,k))^2;
        denomenator = real(x' * (ATilde(:,:,k) + sigma2 * Lambda) * x)
        interference = 0
        for j = 1 : K
            if j ~= k
                interference = interference + abs(H(k,:)*W1(:,j))^2;
            end
        end
        denomenator = interference + sigma2
        
    else
        numerator = real(x'*A(:,:,k)*x);
        denomenator = real(x' * (ATilde(:,:,k) + sigma2 * Lambda) * x);
    end
end

% Wb = beamforming4_4bak(H, sigma2, ATilde,Lambda,a,alpha);
% %Wb = beamforming4_4(H, sigma2)
% Wb = Wb/sqrt(real(trace(Wb'*Wb)));
% Rb(1,:) = calculateRates(H,Wb,sigma2);
% 
% Wb = beamformingMax4_4bak(H, sigma2, ATilde,Lambda,a,alpha);
% % Wb = beamformingMax4_4(H, sigma2)
% Wb = Wb/sqrt(real(trace(Wb'*Wb)));
% Rb(2,:) = calculateRates(H,Wb,sigma2);
% 
% if min(Rb(1,:)) > min(Rb(2,:))
%     selectb = 1;
% else
%     selectb = 2;
% end

[x, xi] = solveFeasible(H,sigma2,ATilde,Lambda,a,alpha);
[xMax, xiMax] = solveFeasibleMax(H,sigma2, ATilde, Lambda,a, alpha);


% W =  H'*inv(H*H'+ alpha*eye(N));
% W = W/sqrt(real(trace(W'*W)));
% R_alpha = calculateRates(H,W,sigma2);
[R_alpha,R1] = calculateRalpha(H, sigma2, alpha,UTEST);

% R = calculateR(x, sigma2, A, ATilde,tTilde2,Lambda,a, UTEST);
W =  H'*inv(H*H'+ U*diag(xi)*U');
W = W/sqrt(real(trace(W'*W)));
RG(1,:)= calculateRates(H,W,sigma2);
xis(1,:) = xi;

W =  H'*inv(H*H'+ U*diag(xiMax)*U');
W = W/sqrt(real(trace(W'*W)));
RG(2,:) = calculateRates(H,W,sigma2);
xis(2,:) = xiMax;

if min(RG(1,:)) >= min(RG(2,:))
    select = 1;
else
    select = 2;
end


sumR_alpha(SINR_SET == SINR) = sumR_alpha(SINR_SET == SINR) + sum(R_alpha);
minR_alpha(SINR_SET == SINR) = minR_alpha(SINR_SET == SINR) + min(R_alpha);

% minRMax(SINR_SET == SINR) = minRMax(SINR_SET == SINR) + min(RMax);
% sumRMax(SINR_SET == SINR) = sumRMax(SINR_SET == SINR) + sum(RMax);

sumR(SINR_SET == SINR) = sumR(SINR_SET == SINR) + sum(RG(select,:));
minR(SINR_SET == SINR) = minR(SINR_SET == SINR) + min(RG(select,:));

% sumRb(SINR_SET == SINR) = sumRb(SINR_SET == SINR) + sum(Rb(select,:));
% minRb(SINR_SET == SINR) = minRb(SINR_SET == SINR) + min(Rb(select,:));
% sumC(SINR_SET == SINR) = sumC(SINR_SET == SINR) + sumCapacity(H, sigma2);

sumAlpha(SINR_SET == SINR) = sumAlpha(SINR_SET == SINR) + K * alpha;
sumXi(SINR_SET == SINR) = sumXi(SINR_SET == SINR) + sum(xis(select,:));
end
end

sumR = sumR/test_num;
minR = minR/test_num;

sumR_alpha = sumR_alpha/test_num;
minR_alpha = minR_alpha/test_num;

% sumRb = sumRb/test_num;
% minRb = minRb/test_num; 
sumAlpha = sumAlpha/test_num;
sumXi = sumXi/test_num;

figure (3)
hold on

plot(SINR_SET, sumRb,'y');
plot(SINR_SET, sumR,'b');
plot(SINR_SET, sumR_alpha, 'r');

figure (2)
hold on
plot(SINR_SET, minRb,'y');
plot(SINR_SET, minR,'b');
plot(SINR_SET, minR_alpha, 'r');


figure (1)
hold on
plot(SINR_SET, sumAlpha,'y');
plot(SINR_SET, sumXi,'b');


save GRCI16_16_1000







    




