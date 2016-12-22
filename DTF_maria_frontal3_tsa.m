function [DTFa]=DTF_maria_frontal3_tsa(low_freq, high_freq, p, Fs, data,thismoment, XYZ, B, name, stemp)
% The data must be in the format nchan x timevector x num_epochs % removes
% ensemble averages - normalization 14-12-2011
nchan=size(data, 1);
num_epochs=size(data,3);

% DTF calculation 
%matlabpool(4)

% Normalization by ensemble mean
for kk=1:nchan
    temp=squeeze(data(kk,:,:));
    ens_mean=mean(temp,2);
    temp_std=std(temp);
    for hh=1:num_epochs
        norm_temp(kk,hh, :)=(temp(:,hh)-ens_mean)./temp_std(hh); % norm_temp is nchan x num_epochs x timeVec
    end
end
% replace data with normalized data The norm_temp einai nchan x num_epoch x
% timevector
data=norm_temp;
% Define frequencies
freq_step=0.5;
fq=low_freq:freq_step:high_freq;
disp(fq(1)); disp('to'); disp(fq(end));
%fq=8:0.25:12;
ORDER_PDC=p;

%*******************************************
% initialization
%*******************************************
avgDC=zeros(nchan,nchan);
avgPDC=zeros(nchan,nchan,length(fq));
avgDTF=zeros(nchan,nchan,length(fq));
ARF=zeros(nchan,ORDER_PDC*nchan);
RCF=zeros(nchan,ORDER_PDC*nchan);
PE=zeros(nchan,ORDER_PDC*nchan+nchan);
DC=zeros(nchan,nchan);

DCtr=zeros(nchan,nchan);
DTFtr=zeros(nchan,nchan,length(fq));
PDCtr=zeros(nchan,nchan,length(fq));
COHtr=zeros(nchan,nchan,length(fq));
for trial=1:num_epochs
    EPOCH=data(1:nchan,trial,:);
    EPOCH=squeeze(EPOCH);map=EPOCH';          
%     for k=1:nrchan
%         %         figure;plot(map(:,k))
%        map(:,k)=map(:,k)-sgolayfilt(map(:,k),3,201);                 
%         %        hold on;plot(map(:,k),'r')
%         %        input('ss')
%         if (ans_normalized=='y')
%               map(:,k)=(map(:,k)-mean(map(:,k)))/std(map(:,k));     
%         else
%               map(:,k)=(map(:,k)-mean(map(:,k)));
%         end
%     end    
   [ARF1,RCF1,PE1,DC1] = mvar_trial(map,ORDER_PDC);
%    [S,h,PDC,COH,DTF,DC,pCOH,dDTF,ffDTF, pCOH2, PDCF]=pdcL(ARF1,RCF1,PE1,ORDER_PDC,fq,Fs);  
%    for k1=1:nrchan
%     for k2=1:nrchan
%         if k1==k2
%         DC(k1,k2)=0;
%         DTF(k1,k2,:)=0;
%         PDC(k1,k2,:)=0;
%         COH(k1,k2,:)=0;
%         end
%     end
%    end
%     DCtr=DCtr+DC;
%     DTFtr=DTFtr+DTF;
%     PDCtr=PDCtr+PDC;
%     COHtr=COHtr+COH; 
    
    % WHAT WE ARE DOING
    ARF=ARF+ARF1;
    RCF=RCF+RCF1;
    PE=PE+ PE1;
    DC=DC+DC1;
    %trial
end
ARF=ARF./trial;
RCF=RCF./trial;
PE=PE./trial;
DC=DC./trial;
[S,h,PDC,COH,DTF,DC,pCOH,dDTF,ffDTF, pCOH2, PDCF]=pdcL(ARF,RCF,PE,ORDER_PDC,fq,Fs);  
% ZERO ON DIAGONAL
for k1=1:nchan
    for k2=1:nchan
        if k1==k2
        DC(k1,k2)=0;
        DTF(k1,k2,:)=0;
        PDC(k1,k2,:)=0;
        COH(k1,k2,:)=0;
        end
    end
end
%figure;
for k1=1:nchan
        for k2=1:nchan
            if k1~=k2
         z=squeeze(abs(pCOH2(k1,k2,:)));
            pCOH_tot(k1,k2)=max(z);
            
             z=squeeze(abs(COH(k1,k2,:)));
            COH_tot(k1,k2)=max(z);
            
%             COH_tot_trial(k1,k2)=max((squeeze(abs(COHtr(k1,k2,:)))));
            z=squeeze(abs(DTF(k1,k2,:)).^2);
           % z(z<0.1)=0;
            DTF_tot(k1,k2)=max(z);
            else
                COH_tot(k1,k2)=0;
                DTF_tot(k1,k2)=0;
            end
        end
end

% SIGNIFICANT VALUES
% PDC_tot(PDC_tot<0.1)=0;
% COH_tot(COH_tot<0.1)=0;

stemp_dir=[thismoment 'DTF-order' num2str(p) '-freqs-' num2str(low_freq) 'to' num2str(high_freq)];
% for jj=1:length(stemp_dir); if (stemp_dir(jj)=='.' || stemp_dir(jj)==' '); stemp_dir(jj)='-';end; end
mkdir(stemp_dir);
cd(stemp_dir);

figure; imagesc(DTF_tot);axis xy; 
set(gca, 'YTickLabel', B); set(gca, 'XTickLabel', B);
set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan);
colorbar('location','EastOutside');
title([name 'DTF-' num2str(low_freq) '-' num2str(high_freq) '-Hz']); 
figure_temp=[name(1:end-4) '-DTF-' num2str(low_freq) '-' num2str(high_freq) '-Hz'];
for jj=1:length(figure_temp); if (figure_temp(jj)=='.' || figure_temp(jj)==':'); figure_temp(jj)='-';end; end
saveas(gcf, figure_temp, 'fig');


% Lines plot
figure; plot2deeg3_frontal(DTF_tot, XYZ, B); title([name(1:end-4) '-' 'DTF-freq-' num2str(low_freq) '-' num2str(high_freq) '-Hz']);
figure_temp=[name(1:3) 'Lines' num2str(low_freq) '-' num2str(high_freq) '-Hz'];
for jj=1:length(figure_temp); if (figure_temp(jj)=='.' || figure_temp(jj)==':'); figure_temp(jj)='-';end; end
saveas(gcf, figure_temp, 'fig');
clear figure_temp 

% Show driving and receiving. Marias
DTF_dr=sum(DTF_tot,1)/29; % DTF driving is the sum on the columns, so sum(c,1), DTF_dr : 1 x29
figure; plot(DTF_dr);set(gca,'Xtick', 1:nchan); set(gca, 'XTickLabel', B);
title([name(1:end-4) '-Driving DTF-freq-' num2str(low_freq) '-' num2str(high_freq) '-Hz']);
figure_temp=[name(1:end-4) '-DTF-driving' num2str(low_freq) '-' num2str(high_freq) '-Hz'];
for jj=1:length(figure_temp); if (figure_temp(jj)=='.' || figure_temp(jj)==':'); figure_temp(jj)='-';end; end
saveas(gcf, figure_temp, 'fig');

DTF_ac=sum(DTF_tot,2)/29;% DTF accepting is the sum on the rows, so sum(c,2), DTF_ac : 1 x29
figure; plot(DTF_ac); set(gca,'Xtick', 1:nchan); set(gca, 'XTickLabel', B);
title([name(1:end-4) '-receiving DTF-freq-' num2str(low_freq) '-' num2str(high_freq) '-Hz']);
figure_temp=[name(1:end-4) '-DTF-receiving' num2str(low_freq) '-' num2str(high_freq) '-Hz'];
for jj=1:length(figure_temp); if (figure_temp(jj)=='.' || figure_temp(jj)==':'); figure_temp(jj)='-';end; end
saveas(gcf, figure_temp, 'fig');
disp([stemp '_figures from DTF_maria_frontal3_tsa.m saved at ' pwd])
% figure;
% subplot(2,2,1);imagesc(COH_tot);colorbar;title('maxCOH')
% subplot(2,2,2);imagesc(DTF_tot);colorbar;title('maxPDC')
% subplot(2,2,3);plot(mean(COH_tot));title('Avg COH');axis tight
% subplot(2,2,4);plot(mean(DTF_tot));title('Driving DTF');axis tight

DTFa.gamma2_all=DTF; % nchan x nchan x numfreqs
DTFa.meangamma=DTF_tot; % nchan x nchan 
DTFa.p=p;
DTFa.low_freq=low_freq;
DTFa.high_freq=high_freq;

 % Can we make it to plot driving and receiving?
% subplot(3,3,1);imagesc(COH_tot_trial);colorbar;title('maxCOH')
% subplot(3,3,2);imagesc(DTF_tot_trial);;colorbar;title('maxDTF')
% subplot(3,3,3);imagesc(PDC_tot_trial);colorbar;title('maxPDC')
% subplot(3,3,4);plot(mean(COH_tot_trial));title('Respose COH');axis tight
% subplot(3,3,5);plot(mean(DTF_tot_trial));title('Respose DTF');axis tight
% subplot(3,3,6);plot(mean(PDC_tot_trial));title('Respose PDC');axis tight
% subplot(3,3,7);plot(mean(COH_tot_trial'));title('Driving COH');axis tight
% subplot(3,3,8);plot(mean(DTF_tot_trial'));title('DRIVING DTF');axis tight
% subplot(3,3,9);plot(mean(PDC_tot_trial'));title('Driving PDC');axis tight

%%
% tic %start old DTF 
% data=norm_temp;
% for jk=1:num_epochs
%        ts=squeeze(data(:,:,jk))'; % timevector x nchan
% %        for cc=1:nchan
% %            ts1(:,cc)=(ts(:,cc)-mean(ts(:,cc)))./std(ts(:,cc));
% %        end
%     gamma2 = DTF(ts,low_freq,high_freq,p,Fs);  %% 10 x 10 x 4
%     %plot(SBC,'DisplayName','SBC','YDataSource','SBC');figure(gcf); title(['chan ' B(ik)])
%     gamma2_all(jk,:,:,:)=gamma2;
%     gamma2=[];
% end
% toc/60
% clear jk ts gamma2 

% for kk=1:nchan,
%        gamma2_all(:,kk,kk,:)=0; %% set zero the diagonal elements 
%      end
% clear kk    
%         
% % gamma = number of epochs x nchan x nchan x 4 (sink, source, frequency index)
% cd(['D:\RIKSHOSPITALET\CNV RIKS\ANALYZED DATASETS\' stemp])
% stemp_dir=[thismoment 'DTF-order' num2str(p) '-freqs-' num2str(low_freq) 'to' num2str(high_freq)];
% mkdir(stemp_dir)
% cd(stemp_dir)
% 
% % make averages for all epochs 
% meanA=squeeze(mean(gamma2_all,1)); % nchan x nchan x Nfreq
% meanFreq=mean(meanA, 3); %% 10 x 10
% gammaFreq=mean(gamma2_all,4); % num_epochs x nchan x nchan 
% meangamma=mean(gammaFreq,1);
% % Matrix DTF plots 
% figure;imagesc(meanFreq);axis xy; 
% set(gca, 'YTickLabel', B); set(gca, 'XTickLabel', B);
% set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan);
% colorbar('location','EastOutside');
% title([name 'DTF: ' num2str(low_freq) '-' num2str(high_freq) ' Hz']); 
% figure_temp=[name(1:end-4) '-matrix-DTF ' num2str(low_freq) '-' num2str(high_freq) '-Hz'];
% saveas(gcf, figure_temp, 'fig')
% clear source
% color_lim=caxis;
% clear figure_temp
% %%\ New Place for Assymetry index
% timeVec2=(1:size(data,2))./Fs;
% 
% % Lines plot
% figure; plot2deeg3_frontal(meanFreq, XYZ, B); title([name '-' 'DTF freq: ' num2str(low_freq) '-' num2str(high_freq) '-Hz']);
% figure_temp=[name(1:end-4) '-Lines-DTF ' num2str(low_freq) '-' num2str(high_freq) '-Hz'];
% saveas(gcf, figure_temp, 'fig')
% clear figure_temp 
% 
% % save all variables into the allready existing DTF substructure
% DTFa.gamma2_all=gamma2_all;
% DTFa.gamma_Freq=gammaFreq;
% DTFa.gamma_time=meanA;
% DTFa.meangamma=squeeze(meangamma);
% DTFa.p=p;
% DTFa.low_freq=low_freq;
% DTFa.high_freq=high_freq;
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

function [CC,NN] = covm(X,Y,Mode);
% COVM generates covariance matrix
% X and Y can contain missing values encoded with NaN.
% NaN's are skipped, NaN do not result in a NaN output. 
% The output gives NaN only if there are insufficient input data
%
% COVM(X,Mode);
%      calculates the (auto-)correlation matrix of X
% COVM(X,Y,Mode);
%      calculates the crosscorrelation between X and Y
%
% Mode = 'M' minimum or standard mode [default]
% 	C = X'*X; or X'*Y correlation matrix
%
% Mode = 'E' extended mode
% 	C = [1 X]'*[1 X]; % l is a matching column of 1's
% 	C is additive, i.e. it can be applied to subsequent blocks and summed up afterwards
% 	the mean (or sum) is stored on the 1st row and column of C
%
% Mode = 'D' or 'D0' detrended mode
%	the mean of X (and Y) is removed. If combined with extended mode (Mode='DE'), 
% 	the mean (or sum) is stored in the 1st row and column of C. 
% 	The default scaling is factor (N-1). 
% Mode = 'D1' is the same as 'D' but uses N for scaling. 
%
% C = covm(...); 
% 	C is the scaled by N in Mode M and by (N-1) in mode D.
% [C,N] = covm(...);
%	C is not scaled, provides the scaling factor N  
%	C./N gives the scaled version. 
%
% see also: DECOVM, XCOVF

%	$Id: covm.m,v 1.21 2005/06/02 15:06:30 schloegl Exp $
%	Copyright (C) 2000-2005 by Alois Schloegl <a.schloegl@ieee.org>	
%       This function is part of the NaN-toolbox
%       http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/NaN/

%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


if nargin<3,
        if nargin==2,
		if isnumeric(Y),
                        Mode='M';
                else
		        Mode=Y;
                        Y=[];
                end;
        elseif nargin==1,
                Mode = 'M';
                Y = [];
        elseif nargin==0,
                error('Missing argument(s)');
        end;
end;        

Mode = upper(Mode);

[r1,c1]=size(X);
if ~isempty(Y)
        [r2,c2]=size(Y);
        if r1~=r2,
                error('X and Y must have the same number of observations (rows).');
                return;
        end;
else
        [r2,c2]=size(X);
end;

if (c1>r1) | (c2>r2),
        warning('Covariance is ill-defined, because of too few observations (rows)');
end;

if ~isempty(Y),
        if (~any(Mode=='D') & ~any(Mode=='E')), % if Mode == M
        	NN = real(X==X)'*real(Y==Y);
	        X(X~=X) = 0; % skip NaN's
	        Y(Y~=Y) = 0; % skip NaN's
        	CC = X'*Y;
        else  % if any(Mode=='D') | any(Mode=='E'), 
	        [S1,N1] = sumskipnan(X,1);
                [S2,N2] = sumskipnan(Y,1);
                
                NN = real(X==X)'*real(Y==Y);
        
	        if any(Mode=='D'), % detrending mode
        		X  = X - ones(r1,1)*(S1./N1);
                        Y  = Y - ones(r1,1)*(S2./N2);
                        if any(Mode=='1'),  %  'D1'
                                NN = NN;
                        else   %  'D0'       
                                NN = max(NN-1,0);
                        end;
                end;
                
                X(X~=X) = 0; % skip NaN's
        	Y(Y~=Y) = 0; % skip NaN's
                CC = X'*Y;
                
                if any(Mode=='E'), % extended mode
                        NN = [r1, N2; N1', NN];
                        CC = [r1, S2; S1', CC];
                end;
	end;
        
else        
        if (~any(Mode=='D') & ~any(Mode=='E')), % if Mode == M
        	tmp = real(X==X);
		NN  = tmp'*tmp;
		X(X~=X) = 0; % skip NaN's
	        CC = X'*X;
        else  % if any(Mode=='D') | any(Mode=='E'), 
	        [S,N] = sumskipnan(X,1);
        	tmp = real(X==X);
                NN  = tmp'*tmp;
                if any(Mode=='D'), % detrending mode
	                X  = X - ones(r1,1)*(S./N);
                        if any(Mode=='1'),  %  'D1'
                                NN = NN;
                        else  %  'D0'      
                                NN = max(NN-1,0);
                        end;
                end;
                
                X(X~=X) = 0; % skip NaN's
                CC = X'*X;
                
                if any(Mode=='E'), % extended mode
                        NN = [r1, N; N', NN];
                        CC = [r1, S; S', CC];
                end;
	end
end;

if nargout<2
        CC = CC./NN; % unbiased
end;
return; 

function [ARF,RCF,PE,DC,varargout] = mvar(Y, Pmax, Mode);
% MVAR estimates Multi-Variate AutoRegressive model parameters
% Several estimation algorithms are implemented, all estimators 
% can handle data with missing values encoded as NaNs.  
%
% 	[AR,RC,PE] = mvar(Y, p);
% 	[AR,RC,PE] = mvar(Y, p, Mode);
%
% INPUT:
%  Y	 Multivariate data series 
%  p     Model order
%  Mode	 determines estimation algorithm 
%
% OUTPUT:
%  AR    multivariate autoregressive model parameter
%  RC    reflection coefficients (= -PARCOR coefficients)
%  PE    remaining error variance
%
% All input and output parameters are organized in columns, one column 
% corresponds to the parameters of one channel.
%
% Mode determines estimation algorithm. 
%  1:  Correlation Function Estimation method using biased correlation function estimation method
%   		also called the "multichannel Yule-Walker" [1,2] 
%  6:  Correlation Function Estimation method using unbiased correlation function estimation method
%
%  2:  Partial Correlation Estimation: Vieira-Morf [2] using unbiased covariance estimates.
%               In [1] this mode was used and (incorrectly) denominated as Nutall-Strand. 
%		Its the DEFAULT mode; according to [1] it provides the most accurate estimates.
%  5:  Partial Correlation Estimation: Vieira-Morf [2] using biased covariance estimates.
%		Yields similar results than Mode=2;
%
%  3:  Partial Correlation Estimation: Nutall-Strand [2] (biased correlation function)
%  7:  Partial Correlation Estimation: Nutall-Strand [2] (unbiased correlation function)
%
% 10:  ARFIT [3] 
% 11:  BURGV [4] 
%
% REFERENCES:
%  [1] A. Schl\"ogl, Comparison of Multivariate Autoregressive Estimators.
%       Signal processing, Elsevier B.V. (in press). 
%       available at http://dx.doi.org/doi:10.1016/j.sigpro.2005.11.007
%  [2] S.L. Marple "Digital Spectral Analysis with Applications" Prentice Hall, 1987.
%  [3] Schneider and Neumaier)
%  [4] Stijn de Waele, 2003
%
%
% A multivariate inverse filter can be realized with 
%   [AR,RC,PE] = mvar(Y,P);
%   e = mvfilter([eye(size(AR,1)),-AR],eye(size(AR,1)),Y);
%  
% see also: MVFILTER, MVFREQZ, COVM, SUMSKIPNAN, ARFIT2

%	$Revision: 1.18 $
%	$Id: mvar.m,v 1.18 2006/04/12 13:11:53 schloegl Exp $
%	Copyright (C) 1996-2006 by Alois Schloegl <a.schloegl@ieee.org>	
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

% CHANGELOG:
% - SUPPORT for ARFIT, BURGV added
% - imputation methods added
% - Correct scaling of PE and covariance matrices. 
% - some experimental versions added. 


% Inititialization
[N,M] = size(Y);

if nargin<2, 
        Pmax=max([N,M])-1;
end;

if iscell(Y)
        Pmax = min(max(N ,M ),Pmax);
        C    = Y;
end;
if nargin<3,
        % according to [1], and other internal validations this is in many cases the best algorithm 
        Mode=2;
end;

[C(:,1:M),n] = covm(Y,'M');
PE(:,1:M)  = C(:,1:M)./n;
if 0,

elseif Mode==6, %%%%% Levinson-Wiggens-Robinson (LWR) algorithm using unbiased correlation function
        C(:,1:M) = C(:,1:M)/N;
        PEF = C(:,1:M);  %?% PEF = PE(:,1:M);
        PEB = C(:,1:M);
        
        for K = 1:Pmax,
                [C(:,K*M+(1:M)),n] = covm(Y(K+1:N,:),Y(1:N-K,:),'M');
                C(:,K*M+(1:M)) = C(:,K*M+(1:M))./n;
		%C{K+1} = C{K+1}/N;

                D = C(:,K*M+(1:M));
                for L = 1:K-1,
                        D = D - ARF(:,L*M+(1-M:0))*C(:,(K-L)*M+(1:M));
                end;
                ARF(:,K*M+(1-M:0)) = D / PEB;	
                ARB(:,K*M+(1-M:0)) = D'/ PEF;	
                for L = 1:K-1,
                        tmp      = ARF(:,L*M+(1-M:0))   - ARF(:,K*M+(1-M:0))*ARB(:,(K-L)*M+(1-M:0));
                        ARB(:,(K-L)*M+(1-M:0)) = ARB(:,(K-L)*M+(1-M:0)) - ARB(:,K*M+(1-M:0))*ARF(:,L*M+(1-M:0));
                        ARF(:,L*M+(1-M:0))   = tmp;
                end;
                
                RCF(:,K*M+(1-M:0)) = ARF(:,K*M+(1-M:0));
                RCB(:,K*M+(1-M:0)) = ARB(:,K*M+(1-M:0));
                
                PEF = [eye(M) - ARF(:,K*M+(1-M:0))*ARB(:,K*M+(1-M:0))]*PEF;
                PEB = [eye(M) - ARB(:,K*M+(1-M:0))*ARF(:,K*M+(1-M:0))]*PEB;
                PE(:,K*M+(1:M)) = PEF;        
        end;
        

elseif Mode==2 %%%%% Partial Correlation Estimation: Vieira-Morf Method [2] with unbiased covariance estimation
	%===== In [1] this algorithm is denoted "Nutall-Strand with unbiased covariance" =====%
        %C(:,1:M) = C(:,1:M)/N;
        F = Y;
        B = Y;
        %PEF = C(:,1:M);
        %PEB = C(:,1:M);
        PEF = PE(:,1:M);
        PEB = PE(:,1:M);
        for K = 1:Pmax,
                [D,n]	= covm(F(K+1:N,:),B(1:N-K,:),'M');
                D = D./n;

		ARF(:,K*M+(1-M:0)) = D / PEB;	
                ARB(:,K*M+(1-M:0)) = D'/ PEF;	
                
                tmp        = F(K+1:N,:) - B(1:N-K,:)*ARF(:,K*M+(1-M:0)).';
                B(1:N-K,:) = B(1:N-K,:) - F(K+1:N,:)*ARB(:,K*M+(1-M:0)).';
                F(K+1:N,:) = tmp;
                
                for L = 1:K-1,
                        tmp      = ARF(:,L*M+(1-M:0))   - ARF(:,K*M+(1-M:0))*ARB(:,(K-L)*M+(1-M:0));
                        ARB(:,(K-L)*M+(1-M:0)) = ARB(:,(K-L)*M+(1-M:0)) - ARB(:,K*M+(1-M:0))*ARF(:,L*M+(1-M:0));
                        ARF(:,L*M+(1-M:0))   = tmp;
                end;
                
                RCF(:,K*M+(1-M:0)) = ARF(:,K*M+(1-M:0));
                RCB(:,K*M+(1-M:0)) = ARB(:,K*M+(1-M:0));
                
                [PEF,n] = covm(F(K+1:N,:),F(K+1:N,:),'M');
                PEF = PEF./n;

		[PEB,n] = covm(B(1:N-K,:),B(1:N-K,:),'M');
                PEB = PEB./n;

                PE(:,K*M+(1:M)) = PEF;        
        end;
        



end;


if any(ARF(:)==inf),
% Test for matrix division bug. 
% This bug was observed in LNX86-ML5.3, 6.1 and 6.5, but not in SGI-ML6.5, LNX86-ML6.5, Octave 2.1.35-40; Other platforms unknown.
p = 3;
FLAG_MATRIX_DIVISION_ERROR = ~all(all(isnan(repmat(0,p)/repmat(0,p)))) | ~all(all(isnan(repmat(nan,p)/repmat(nan,p))));

if FLAG_MATRIX_DIVISION_ERROR, 
	%fprintf(2,'### Warning MVAR: Bug in Matrix-Division 0/0 and NaN/NaN yields INF instead of NAN.  Workaround is applied.\n');
	warning('MVAR: bug in Matrix-Division 0/0 and NaN/NaN yields INF instead of NAN.  Workaround is applied.');

	%%%%% Workaround 
	ARF(ARF==inf)=NaN;
	RCF(RCF==inf)=NaN;
end;
end;	

%MAR   = zeros(M,M*Pmax);
DC     = zeros(M);
for K  = 1:Pmax,
%       VAR{K+1} = -ARF(:,K*M+(1-M:0))';
        DC  = DC + ARF(:,K*M+(1-M:0))'.^2; %DC meausure [3]
end;
% figure;imagesc(DC)
% input('ss')
% Developers:  What's New (Koders Blog) | Getting Started | Add Projects | Downloads
% Add Search to Your Site | Open Source Zeitgeist | Advertise | About Koders
% 
% RSS Feed
% 
% Copyright ? 2006 Koders  -  Searching 512,549,710 lines of code.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [ARF,RCF,PE,DC,varargout] = mvar(Y, Pmax, Mode);
function [ARF,RCF,PE,DC,varargout] = mvar_trial(Y, Pmax, Mode);
% Inititialization
[N,M] = size(Y);
if iscell(Y)
        Pmax = min(max(N ,M ),Pmax);
        C    = Y;
end;
if nargin<3,
        % according to [1], and other internal validations this is in many cases the best algorithm 
        Mode=2;
end;
[C(:,1:M),n] = covm(Y,'M');
PE(:,1:M)  = C(:,1:M)./n;
if 0,
elseif Mode==2 %%%%% Partial Correlation Estimation: Vieira-Morf Method [2] with unbiased covariance estimation
	%===== In [1] this algorithm is denoted "Nutall-Strand with unbiased covariance" =====%
        %C(:,1:M) = C(:,1:M)/N;
        F = Y;
        B = Y;
        %PEF = C(:,1:M);
        %PEB = C(:,1:M);
        PEF = PE(:,1:M);
        PEB = PE(:,1:M);
        for K = 1:Pmax,
                [D,n]	= covm(F(K+1:N,:),B(1:N-K,:),'M');
                D = D./n;

		       ARF(:,K*M+(1-M:0)) = D / PEB;	
                ARB(:,K*M+(1-M:0)) = D'/ PEF;	
                
                tmp        = F(K+1:N,:) - B(1:N-K,:)*ARF(:,K*M+(1-M:0)).';
                B(1:N-K,:) = B(1:N-K,:) - F(K+1:N,:)*ARB(:,K*M+(1-M:0)).';
                F(K+1:N,:) = tmp;
                
                for L = 1:K-1,
                        tmp      = ARF(:,L*M+(1-M:0))   - ARF(:,K*M+(1-M:0))*ARB(:,(K-L)*M+(1-M:0));
                        ARB(:,(K-L)*M+(1-M:0)) = ARB(:,(K-L)*M+(1-M:0)) - ARB(:,K*M+(1-M:0))*ARF(:,L*M+(1-M:0));
                        ARF(:,L*M+(1-M:0))   = tmp;
                end;
                
                RCF(:,K*M+(1-M:0)) = ARF(:,K*M+(1-M:0));
                RCB(:,K*M+(1-M:0)) = ARB(:,K*M+(1-M:0));
                
                [PEF,n] = covm(F(K+1:N,:),F(K+1:N,:),'M');
                PEF = PEF./n;

		[PEB,n] = covm(B(1:N-K,:),B(1:N-K,:),'M');
                PEB = PEB./n;

                PE(:,K*M+(1:M)) = PEF;        
        end;
        



end;


if any(ARF(:)==inf),
% Test for matrix division bug. 
% This bug was observed in LNX86-ML5.3, 6.1 and 6.5, but not in SGI-ML6.5, LNX86-ML6.5, Octave 2.1.35-40; Other platforms unknown.
p = 3;
FLAG_MATRIX_DIVISION_ERROR = ~all(all(isnan(repmat(0,p)/repmat(0,p)))) | ~all(all(isnan(repmat(nan,p)/repmat(nan,p))));

if FLAG_MATRIX_DIVISION_ERROR, 
	%fprintf(2,'### Warning MVAR: Bug in Matrix-Division 0/0 and NaN/NaN yields INF instead of NAN.  Workaround is applied.\n');
	warning('MVAR: bug in Matrix-Division 0/0 and NaN/NaN yields INF instead of NAN.  Workaround is applied.');

	%%%%% Workaround 
	ARF(ARF==inf)=NaN;
	RCF(RCF==inf)=NaN;
end;
end;	

%MAR   = zeros(M,M*Pmax);
DC     = zeros(M);
for K  = 1:Pmax,
%       VAR{K+1} = -ARF(:,K*M+(1-M:0))';
        DC  = DC + ARF(:,K*M+(1-M:0))'.^2; %DC meausure [3]
end;
% figure;imagesc(DC)
% input('ss')
% Developers:  What's New (Koders Blog) | Getting Started | Add Projects | Downloads
% Add Search to Your Site | Open Source Zeitgeist | Advertise | About Koders
% 
% RSS Feed
% 
% Copyright ? 2006 Koders  -  Searching 512,549,710 lines of code.
function data=filter_data(data,from,to,Fs)
% eg from=5; to=45
%Fs=1000;
[xr,xc]=size(data);
nrchan=min(xr,xc);
t=1:length(data(1,:))';
% figure;subplot(1,2,1);
% plot(data);axis tight;
% [p,f]=pwelch(data(1,:),[],[],[],Fs);
% subplot(1,2,2);hold on;plot(f,10*log10(p))
% axis tight


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% band pass filtering data in 5 -45
for k=1:nrchan
    data(k,:) = eegfiltfft(data(k,:), Fs, from, to);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(1,2,1);hold on;plot(data(1,:),'r');legend('orig','filtered')
% axis tight;xlabel('time');ylabel('Ampl');title('Chan 1')
% 
% % [p,f]=pwelch(data(1,:),[],[],[],Fs);
% % subplot(1,2,2);hold on;plot(f,10*log10(p),'r');axis tight
% axis tight; xlabel('frequency [Hz]');ylabel('PSD')

% subplot(1,2,1);imagesc(x,y,Zi);axis tight
% subplot(1,2,2);plot(y, x, 'x', 'Color', 'black', 'markersize', 5); hold on
% contour(Xi, Yi, sumLaplac2D); hold off;axis tight
% input('')
%return;


