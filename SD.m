function SD(Y, C, N)

global bEst;
global radius;
bEst = 0;
[Q,R] = qr(C);
a_hat_SD_ZF = sign(C\Y);
radius = (norm(Y - C*a_hat_SD_ZF).^2)*1.2;
z = Q' * Y;
b_hat = zeros(N,1);
funct(z, R, [1 -1], N , 0, b_hat);

clear bEst;