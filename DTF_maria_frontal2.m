function [DTFa]=DTF_maria_frontal2(low_freq, high_freq, p, Fs, data,thismoment, XYZ, B, name, stemp)
% The data must be in the format nchan x timevector x num_epochs % removes
% ensemble averages - normalization 14-12-2011
nchan=size(data, 1);
num_epochs=size(data,3);

% DTF calculation 
%matlabpool(4)

% normalization by ensemble mean
for kk=1:nchan
    temp=squeeze(data(kk,:,:));
    ens_mean=mean(temp,2);
    temp_std=std(temp);
    for hh=1:num_epochs
        norm_temp(kk,:, hh)=(temp(:,hh)-ens_mean)./temp_std(hh);
    end
    
end
tic
data=norm_temp;
for jk=1:num_epochs
       ts=squeeze(data(:,:,jk))'; % timevector x nchan
%        for cc=1:nchan
%            ts1(:,cc)=(ts(:,cc)-mean(ts(:,cc)))./std(ts(:,cc));
%        end
    gamma2 = DTF(ts,low_freq,high_freq,p,Fs);  %% 10 x 10 x 4
    %plot(SBC,'DisplayName','SBC','YDataSource','SBC');figure(gcf); title(['chan ' B(ik)])
    gamma2_all(jk,:,:,:)=gamma2;
    gamma2=[];
end
toc/60
clear jk ts gamma2 

for kk=1:nchan,
       gamma2_all(:,kk,kk,:)=0; %% set zero the diagonal elements 
     end
clear kk    
        
% gamma = number of epochs x nchan x nchan x 4 (sink, source, frequency index)
cd(['D:\RIKSHOSPITALET\CNV RIKS\ANALYZED DATASETS\' stemp])
stemp_dir=[thismoment 'DTF-order' num2str(p) '-freqs-' num2str(low_freq) 'to' num2str(high_freq)];
mkdir(stemp_dir)
cd(stemp_dir)

% make averages for all epochs 
meanA=squeeze(mean(gamma2_all,1)); % nchan x nchan x Nfreq
meanFreq=mean(meanA, 3); %% 10 x 10
gammaFreq=mean(gamma2_all,4); % num_epochs x nchan x nchan 
meangamma=mean(gammaFreq,1);
% Matrix DTF plots 
figure;imagesc(meanFreq);axis xy; 
set(gca, 'YTickLabel', B); set(gca, 'XTickLabel', B);
set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan);
colorbar('location','EastOutside');
title([name 'DTF: ' num2str(low_freq) '-' num2str(high_freq) ' Hz']); 
figure_temp=[name(1:end-4) '-matrix-DTF ' num2str(low_freq) '-' num2str(high_freq) '-Hz'];
saveas(gcf, figure_temp, 'fig')
clear source
color_lim=caxis;
clear figure_temp
%%\ New Place for Assymetry index
timeVec2=(1:size(data,2))./Fs;

% Lines plot
figure; plot2deeg3_frontal(meanFreq, XYZ, B); title([name '-' 'DTF freq: ' num2str(low_freq) '-' num2str(high_freq) '-Hz']);
figure_temp=[name(1:end-4) '-Lines-DTF ' num2str(low_freq) '-' num2str(high_freq) '-Hz'];
saveas(gcf, figure_temp, 'fig')
clear figure_temp 

% save all variables into the allready existing DTF substructure
DTFa.gamma2_all=gamma2_all;
DTFa.gamma_Freq=gammaFreq;
DTFa.gamma_time=meanA;
DTFa.meangamma=squeeze(meangamma);
DTFa.p=p;
DTFa.low_freq=low_freq;
DTFa.high_freq=high_freq;