function rate = RATE(X,N,Es)

NFFT = 2^nextpow2(N);
Y = fft(X,NFFT)/N;

qu = 2*abs(Y(1:NFFT/2+1));
% %% XCIN
% ANCGP = qu(42);    
% ANC2P = qu(83); 
% ANC3P = qu(222);    
% ANC4P = qu(443);    
% ANC5P = qu(664);    
% t = [ANCGP ANC2P ANC3P ANC4P ANC5P];

%% XOUTER
%ANCGP = qu(7);                              
ANC2P = qu(22);    
ANC3P = qu(43);    
ANC4P = qu(64); 
ANC5P = qu(84);
t = [ANC2P ANC3P ANC4P ANC5P];

rate =  sum(t.^2)/ sum(qu(2:1:Es).^2);
end