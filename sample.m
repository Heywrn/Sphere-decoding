clc
clear

N = 50;
[G] = FTN (N);    % G = toeplitz matrix
load('matrixR');  % R = randn, Hermition matrix


startSNR = 0;     % start of SNR(dB)
stepSNR = 2;      % increment of SNR
endSNR = 18;      % end of SNR (dB)
SNR=[startSNR:stepSNR:endSNR];  % SNR range
Eb_N0_dB  = [startSNR:stepSNR:endSNR]; % multiple Es/N0 values
%Es_N0_dB  = Eb_N0_dB + 10*log10(km);
BER_FSD_G=zeros(N,length(SNR));

N_sim = 10;
for m=1:length(SNR)
    m
    %a=sqrt(10^(SNR(m)/10))/sqrt(L);
    noiseV = 10^(-Eb_N0_dB(m)/20);
    
    for mm=1:N_sim


   
B=sign(randn(N,1));
Noise = randn(N,1);

YG = G*B + noiseV*Noise;
        

%% FSD
        
GG = [ -G'*G, G'*YG; YG'*G,1];
        
[Vg, Dg, Ug] = eig(GG);
GGN = GG+ 2*max(diag(Dg))*eye(N+1);
AGG = chol(GGN);


y = YG;
A = AGG;
BFSD = sign(randn(N,1));

s = zeros(N+1);
% Initialize s;
L = N;
s(1,L) = 1;

j1 = 1;

% Fisrt, get the s(1)
while(j1<L)
    mid = 0;
    for i=L:-1:L-j1+1
        mid = mid + A(:,L-j1)' * A(:,i) * s(1,i);
    end
    s(1,L-j1) = sign(mid);
    j1 = j1+1;
end


% Based on s(1), we can get s(n)
n = 2;

while(n<=L)
    s(n,L-n+1) = -s(1,L-n+1);
    for i=L-n+2:L
        s(n,i) = s(1,i);
    end
    j = n;
    while(j<L)
        mid = 0;
        for i=L:-1:L-j+1
            mid = mid + A(:,L-j)' * A(:,i) * s(n,i);
        end
        s(n,L-j) = sign(mid);
        j = j+1;
    end
    n = n+1;
end

% calculate which s is the best one



be1 = norm(y - A(1:50,1:50) * s(1:50,1))^2;
be2 = norm(y + A(1:50,1:50) * s(1:50,1))^2;
if(be1 > be2)
    be = be1;
else
    be = be2;
end
BFSD = s(1:50,1);

for q = 2:L
    Cal1 = norm(y - A(1:50,1:50) * s(1:50,q))^2;
    Cal2 = norm(y + A(1:50,1:50) * s(1:50,q))^2;
    if(Cal1 > Cal2)
        Cal = Cal1;
    else
        Cal = Cal2;
    end
    if(Cal > be)
        be = Cal;
        BFSD = s(1:50,q);
    end
end


BER_FSD_G(:,m) =BER_FSD_G(:,m) + (abs(BFSD -B));
    end
end
BER_FSD_G_max=max(BER_FSD_G)/N_sim ;
BER_FSD_G_ave=sum(BER_FSD_G)/N_sim/N;



