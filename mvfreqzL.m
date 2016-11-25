function [S,h,PDC,COH,DTF,DC,pCOH,dDTF,ffDTF, pCOH2, PDCF, coh]=mvfreqzL(B,A,C,N,Fs)
% MVFREQZ multivariate frequency response
% [S,h,PDC,COH,DTF,DC,pCOH,dDTF,ffDTF,pCOH2,PDCF] = mvfreqz(B,A,C,N,Fs)
%
% INPUT: 
% ======= 
% A, B	multivariate polynomials defining the transfer function
%
%    a0*Y(n) = b0*X(n) + b1*X(n-1) + ... + bq*X(n-q)
%                          - a1*Y(n-1) - ... - ap*Y(:,n-p)
%
%  A=[a0,a1,a2,...,ap] and B=[b0,b1,b2,...,bq] must be matrices of
%  size  Mx((p+1)*M) and Mx((q+1)*M), respectively. 
%
%  C is the covariance of the input noise X
%  N if scalar, N is the number of frequencies 
%    if N is a vector, N are the designated frequencies. 
%  Fs sampling rate [default 2*pi]
% 
%  A,B and C can by obtained from a multivariate time series 
%       through the following commands: 
%  [AR,RC,PE] = mvar(Y,P);
%       M = size(AR,1); % number of channels       
%       A = [eye(M),-AR];
%       B = eye(M); 
%       C = PE(:,M*P+1:M*(P+1)); 
%
% OUTPUT: 
% ======= 
% S   	power spectrum
% PDC 	partial directed coherence
% DC  	directed coupling	
% COH 	coherency (complex coherence)
% DTF 	directed transfer function
% pCOH 	partial coherence
% dDTF 	direct Directed Transfer function
% ff[DTF full frequency Directed Transfer Function 
% pCOH2 partial coherence -alternative method 
%
%
% see also: FREQZ, MVFILTER, MVAR
%
% 
% REFERENCE(S):
% H. Liang et al. Neurocomputing, 32-33, pp.891-896, 2000. 
% L.A. Baccala and K. Samashima, Biol. Cybern. 84,463-474, 2001. 
% A. Korzeniewska, et al. Journal of Neuroscience Methods, 125, 195-207, 2003. 
% Piotr J. Franaszczuk, Ph.D. and Gregory K. Bergey, M.D.
% 	Fast Algorithm for Computation of Partial Coherences From Vector Autoregressive Model Coefficients
%	World Congress 2000, Chicago. 
% Nolte G, Bai O, Wheaton L, Mari Z, Vorbach S, Hallett M.
%	Identifying true brain interaction from EEG data using the imaginary part of coherency.
%	Clin Neurophysiol. 2004 Oct;115(10):2292-307. 
% Schlgl A., Supp G.
%       Analyzing event-related EEG data with multivariate autoregressive parameters.
%       (Eds.) C. Neuper and W. Klimesch, 
%       Progress in Brain Research: Event-related Dynamics of Brain Oscillations. 
%       Analysis of dynamics of brain oscillations: methodological advances. Elsevier. 


%	$Id: mvfreqz.m,v 1.8 2006/05/04 14:36:37 schloegl Exp $
%	Copyright (C) 1996-2005 by Alois Schloegl <a.schloegl@ieee.org>	
%       This is part of the TSA-toolbox. See also 
%       http://hci.tugraz.at/schloegl/matlab/tsa/
%       http://octave.sourceforge.net/
%       http://biosig.sourceforge.net/


% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.

[K1,K2] = size(A);
p = K2/K1-1;
%a=ones(1,p+1);
[K1,K2] = size(B);
q = K2/K1-1;
%b=ones(1,q+1);
if nargin<3
        C = eye(K1,K1);
end;
if nargin<4,
        N = 512;
end;
if nargin<5,
        Fs= 1;        
end;
if all(size(N)==1),	
        f = (0:N-1)/N;
else
        f = N;
        N = length(N);
end;
s = exp(i*2*pi*f/Fs);
z = i*2*pi/Fs;

h=zeros(K1,K1,N);
g=zeros(K1,K1,N);
S=zeros(K1,K1,N);
S1=zeros(K1,K1,N);
DTF=zeros(K1,K1,N);
COH=zeros(K1,K1,N);
%COH2=zeros(K1,K1,N);
PDC=zeros(K1,K1,N);
PDCF = zeros(K1,K1,N);
pCOH = zeros(K1,K1,N);
invC=pinv(C);
tmp1=zeros(1,K1);
tmp2=zeros(1,K1);

M = zeros(K1,K1,N);
detG = zeros(N,1);

for n=1:N,
        atmp = zeros(K1);
        for k = 1:p+1,
                atmp = atmp + A(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
        end;        
        btmp = zeros(K1);
        for k = 1:q+1,
                btmp = btmp + B(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
        end;        
        h(:,:,n) = atmp\btmp;        
        S(:,:,n) = h(:,:,n)*C*h(:,:,n)';        
        S1(:,:,n) = h(:,:,n)*h(:,:,n)';        
        
        for k1 = 1:K1,
                tmp = squeeze(atmp(:,k1));
                tmp1(k1) = sqrt(tmp'*tmp);
                tmp2(k1) = sqrt(tmp'*invC*tmp);
        end;
        
        PDCF(:,:,n) = abs(atmp)./tmp2(ones(1,K1),:);
        PDC(:,:,n)  = abs(atmp)./tmp1(ones(1,K1),:);
        
        g = atmp/btmp;        
        G(:,:,n) = g'*invC*g;
        detG(n) = det(G(:,:,n));        
        
end;

if nargout<4, return; end;

%%%%% directed transfer function
for k1=1:K1;
        DEN=sum(abs(h(k1,:,:)).^2,2);	        
        for k2=1:K2;
                %COH2(k1,k2,:) = abs(S(k1,k2,:).^2)./(abs(S(k1,k1,:).*S(k2,k2,:)));
                COH(k1,k2,:) = (S(k1,k2,:))./sqrt(abs(S(k1,k1,:).*S(k2,k2,:)));
                coh(k1,k2,:) = (S1(k1,k2,:))./sqrt(abs(S1(k1,k1,:).*S1(k2,k2,:)));
                %DTF(k1,k2,:) = sqrt(abs(h(k1,k2,:).^2))./DEN;	        
%                 DTF(k1,k2,:) = abs(h(k1,k2,:)).^2./sqrt(DEN);
                DTF(k1,k2,:) = abs(h(k1,k2,:)).^2./(DEN);
                ffDTF(k1,k2,:) = abs(h(k1,k2,:))./sqrt(sum(DEN,3));
                pCOH2(k1,k2,:) = abs(G(k1,k2,:).^2)./(G(k1,k1,:).*G(k2,k2,:));
                
                M(k2,k1,:) = ((-1)^(k1+k2))*squeeze(G(k1,k2,:))./detG; % oder ist M = G?
        end;
end;

for k1=1:K1;
        for k2=1:K2;
                pCOH(k1,k2,:) = abs(M(k1,k2,:).^2)./(M(k1,k1,:).*M(k2,k2,:));
        end;
end;

dDTF = pCOH2.*ffDTF; 

if nargout<6, return; end;

DC = zeros(K1);
for k = 1:p,
        DC = DC + A(:,k*K1+(1:K1)).^2;
end;        
%DC=DC./sum(sum(A.^2));

if nargout<7, return; end;


for n=1:N,
        %COH2(k1,k2,:) = abs(S(k1,k2,:).^2)./(abs(S(k1,k1,:).*S(k2,k2,:)));
        M(k1,k2,n) = det(squeeze(S([1:k1-1,k1+1:K1],[1:k2-1,k2+1:K2],n)));
end;

for k1=1:K1;
        for k2=1:K2;
                for n=1:N,
                        %COH2(k1,k2,:) = abs(S(k1,k2,:).^2)./(abs(S(k1,k1,:).*S(k2,k2,:)));
                        M(k1,k2,n) = det(squeeze(S([1:k1-1,k1+1:K1],[1:k2-1,k2+1:K2],n)));
                end;
        end;
end;

for k1=1:K1;
        for k2=1:K2;
                pCOH(k1,k2,:) = abs(M(k1,k2,:).^2)./(M(k1,k1,:).*M(k2,k2,:));
        end;
end;
