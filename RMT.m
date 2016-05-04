% This program simulate the CoMP conference paper.

clear all
DEBUG = 0

t0 = clock;
%Parameter: # of MS, # of BS, # of receiving antenna, # of transmitting
%antenna.
N = 4; 
K = 4;

%%%%%%%%%%%%%%%%%%%Generate h(:,k)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1: K
%     l(k) = rand()
%     generatePathLoss
    l(k) = 1
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
    x = (randn(1,N)+sqrt(-1)*randn(1,N))/sqrt(2);
    h(k,:) = l(k) * x;
end

for rho = 0.0001:0.0001:1
    m = calculateM(l, rho, K, K/N);
    g = calculateg(l, rho, K, K/N);
    dm = calculatedMdrho(l, m, K, K/N);
    dg = calculatedg(l, rho, K, K/N);
end

rho_set = 0.001:0.001:1;
sigma2 = 0.1;
for i = rho_set
    m(rho_set == i)  = calculateM(l, i, K, K/N);
    dm(rho_set == i) = calculatedMdrho(l, i, K, K/N);
    R_inf(rho_set == i) = calculateRinf(l, i, K, K/N,sigma2);
    Ru_inf(rho_set == i) = calculateRinfu(l, i, K, K/N,sigma2);
    R(rho_set == i) = calculateR(h, sigma2, i * N);
end

index = find(R_inf  == max(R_inf));
rho_set(index)
figure
plot(rho_set,m)
figure
plot(rho_set,dm)
figure
plot(rho_set,real(R_inf), '-sb')
hold on
plot(rho_set,real(Ru_inf),'-dr')
plot(rho_set,real(R),'-og')


    




