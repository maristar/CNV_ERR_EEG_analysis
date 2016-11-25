function [S,h,PDC,COH,DTF,DC,pCOH,dDTF,ffDTF, pCOH2, PDCF, coh]=pdcL(AR,RC,PE,P,N,Fs)
% N - nr of frequencies, or fq vector
% Fs sampling fq
% Y- the mltivariate data A matrix of data points x nr channels
% P- order of the mvar model

%[] = mvar(Y,P);
% size(AR)
M = size(AR,1); % number of channels       
A = [eye(M),-AR];
B = eye(M); 
C = PE(:,M*P+1:M*(P+1)); 
[S,h,PDC,COH,DTF,DC,pCOH,dDTF,ffDTF, pCOH2, PDCF, coh]=mvfreqzL(B,A,C,N,Fs);