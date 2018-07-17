function [G] = FTN (N)

alpha = 0.3;
tau = 0.5;

M = 2;
N_back = 3;           % number of symbols to go back and correct for SSD_GBK

% interference
fd = 1;
fs = 100;
gdelay = 4;
h = rcosine(fd,fs,'sqrt',alpha,gdelay); % group delay(samples) = gdelay*fs+1
h_rc = conv(h,h);
pulse_loc  = 2*gdelay*fs + 1;
g = h_rc(pulse_loc : tau*fs : end);
% g_short = g(1:10);
g_short = g;
if length(g_short) > N
    disp('   Warning: increase N     ');
end
ISI = [g_short zeros(1,N-length(g_short))];
for jj = 1:length(ISI)
    II(jj,:) = circshift(ISI',jj-1); % II is the ISI matrix
end
II = II - tril(II,-1);
II = triu(II)+triu(II,1)';
G = II;