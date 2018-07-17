close all
clear all


N = 18;
[G] = FTN (N);    % G = toeplitz matrix
load('matrixR');  % R = randn, Hermition matrix
load('matrixF');  % F = Orthogonal matrix (eig(R))


startSNR = 0;     % start of SNR(dB)
stepSNR = 2;      % increment of SNR
endSNR = 18;      % end of SNR (dB)
SNR=[startSNR:stepSNR:endSNR];  % SNR range
Eb_N0_dB  = [startSNR:stepSNR:endSNR]; % multiple Es/N0 values
%Es_N0_dB  = Eb_N0_dB + 10*log10(km);

% QF Factorization
[q r]=qr(G);
Q1=q(1:N,1:N);
%Q2=q(1:m-n,1:m);
R=r(1:N,1:N);



N_sim = 100;

start = clock


BER_MF_G=zeros(N,length(SNR));
BER_MF_R=zeros(N,length(SNR));
BER_MF_F=zeros(N,length(SNR));

BER_SD_G=zeros(N,length(SNR));
BER_SD_R=zeros(N,length(SNR));
BER_SD_F=zeros(N,length(SNR));

BER_FSD_G=zeros(N,length(SNR));
BER_FSD_R=zeros(N,length(SNR));
BER_FSD_F=zeros(N,length(SNR));


for m=1:length(SNR)
    m
    %a=sqrt(10^(SNR(m)/10))/sqrt(L);
    noiseV = 10^(-Eb_N0_dB(m)/20);
    
    for mm=1:N_sim
        
        
        B=sign(randn(N,1));
        Noise = randn(N,1);
        %noiseV = 0;
        YG = G*B + noiseV*Noise;
        YR = R*B + noiseV*Noise;
        YF = F*B + noiseV*Noise;

        [bGMF] = sign(G'*YG);  % matched filter (MF)
        [bRMF] = sign(R'*YR);  % matched filter (MF)
        [bFMF] = sign(F'*YF);  % matched filter (MF)
        
        %bGSD = SD(G, YG, [1 -1]); % sphered decoder
        %bRSD = SD(R, YR, [1 -1]); % sphered decoder
        %bFSD = SD(F, YF, [1 -1]); % sphered decoder

        % [bGFSD] = FSD(YG, G, N); %Finite Step SD
        %[bRFSD] = FSD(YR, R, N ); %Finite Step SD
        %[bFFSD] = FSD(YF, F, N ); %Finite Step SD
        
        BER_MF_G(:,m) =BER_MF_G(:,m) + (abs(bGMF -B));
        BER_MF_R(:,m) =BER_MF_R(:,m) + (abs(bRMF -B));
        BER_MF_F(:,m) =BER_MF_F(:,m) + (abs(bFMF -B));
        
%         if (length(bGSD)~= 1)
%             BER_SD_G(:,m) =BER_SD_G(:,m) + (abs(bGSD -B));
%         end
%         if (length(bRSD)~= 1)
%             BER_SD_R(:,m) =BER_SD_R(:,m) + (abs(bRSD -B));
%         end
%         if (length(bFSD)~= 1)
%             BER_SD_F(:,m) =BER_SD_F(:,m) + (abs(bFSD -B));
%         end
        %
        %         BER_FSD_G(:,m) =BER_FSD_G(:,m) + (abs(bGFSD -B));
        %         BER_FSD_R(:,m) =BER_FSD_R(:,m) + (abs(bRFSD -B));
        %         BER_FSD_F(:,m) =BER_FSD_F(:,m) + (abs(bFFSD -B));
    end
end


BER_MF_G_max=max(BER_MF_G)/N_sim ;
BER_MF_R_max=max(BER_MF_R)/N_sim ;
BER_MF_F_max=max(BER_MF_F)/N_sim ;

BER_SD_G_max=max(BER_SD_G)/N_sim ;
BER_SD_R_max=max(BER_SD_R)/N_sim ;
BER_SD_F_max=max(BER_SD_F)/N_sim ;

BER_FSD_G_max=max(BER_FSD_G)/N_sim ;
BER_FSD_R_max=max(BER_FSD_R)/N_sim ;
BER_FSD_F_max=max(BER_FSD_F)/N_sim ;

BER_MF_G_ave=sum(BER_MF_G)/N_sim/N;
BER_MF_R_ave=sum(BER_MF_R)/N_sim/N;
BER_MF_F_ave=sum(BER_MF_F)/N_sim/N;

BER_SD_G_ave=sum(BER_SD_G)/N_sim/N;
BER_SD_R_ave=sum(BER_SD_R)/N_sim/N;
BER_SD_F_ave=sum(BER_SD_F)/N_sim/N;

BER_FSD_G_ave=sum(BER_FSD_G)/N_sim/N;
BER_FSD_R_ave=sum(BER_FSD_R)/N_sim/N;
BER_FSD_F_ave=sum(BER_FSD_F)/N_sim/N;


endtime =  clock;
elapsed =  (endtime - start)*[0 0 24*60^2 60.^[2 1 0]]'


figure
semilogy( SNR, BER_MF_G_ave,  '-*r', SNR, BER_SD_G_ave,  '-xr', SNR, BER_FSD_G_ave,  '-or',...
    SNR, BER_MF_R_ave,  '-*b', SNR, BER_SD_R_ave,  '-xb', SNR, BER_FSD_R_ave,  '-ob',...
    SNR, BER_MF_F_ave,  '-*g', SNR, BER_SD_F_ave,  '-xg', SNR, BER_FSD_F_ave,  '-og')
legend( 'MF G',  'SD G', 'FSD G', 'MF R',  'SD R', 'FSD R', 'MF F',  'SD F', 'FSD F')
grid on
ylabel('Average BER')
xlabel('SNR in dB')


figure
semilogy( SNR, BER_MF_G_max,  '-*r', SNR, BER_SD_G_max,  '-xr', SNR, BER_FSD_G_max,  '-or',...
    SNR, BER_MF_R_max,  '-*b', SNR, BER_SD_R_max,  '-xb', SNR, BER_FSD_R_max,  '-ob',...
    SNR, BER_MF_F_max,  '-*g', SNR, BER_SD_F_max,  '-xg', SNR, BER_FSD_F_max,  '-og')
legend( 'MF G',  'SD G', 'FSD G', 'MF R',  'SD R', 'FSD R', 'MF F',  'SD F', 'FSD F')
grid on
ylabel('MAX BER')
xlabel('SNR in dB')








