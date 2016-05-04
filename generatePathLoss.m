% This function is to calculate the Stiljes trasform of -\rho.


function l = generatePathLoss(T)
r = 1;
r_min = 0.001;
c = r * r - r_min * r_min;

u = rand(T,1);
l = (c * u + r^2).^(-1/2);

end