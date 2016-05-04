function [x, xi] = solveFeasible(sigma2,ATilde,Lambda,a)
K = size(Lambda,1);
t_max = 10;
t_min = 0;
threshold = 0.001

while t_max - t_min > threshold

t = (t_max + t_min)/2;
cvx_begin quiet
cvx_precision best
variable x(K)
maximize(sum(x))
for k = 1 : K
    norm(t*( ATilde(:,:,k) + sigma2 * Lambda)^(1/2) *x) <=  real(a(:,k,k))'* x
end
subject to
for k = 1 : K
  x(k) <= Lambda(k,k)^(-1)
  x(k) >= 0.001
end
cvx_end
if strcmp(cvx_status, 'Infeasible') == 1 || strcmp(cvx_status, 'Failed') 
    t_max = t;
else
    t_min = t;
    result = x;
end

for k = 1 : K
    xi(k) = 1/x(k) - Lambda(k,k) ;
end
x= x;
end