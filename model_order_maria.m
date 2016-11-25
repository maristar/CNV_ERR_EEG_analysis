function p=model_order_maria(dataf,pstart, pend)
%% model order for DTF
%% it uses the dataf so the data is in the format nchan x single trial duration x number of single trials
nchan=size(dataf,1);
num_epochs=size(dataf,3);
%load B10
% pstart=1;
% pend=25;
N=pend-pstart+1;
for ik=1:nchan
    for jk=1:num_epochs
    [popt n,m, C,SBC,FPE,th]=arfit_m(dataf(ik,:,jk)',pstart,pend);
    %figure; plot(SBC,'DisplayName','SBC','YDataSource','SBC');figure(gcf); title(['chan ' B(ik)])
    SBC_all(ik,jk,:)=SBC;
    FPE_all(ik,jk,:)=FPE;
    %popt_all(ik,jk,:)=popt;
    end
end

for ik=1:nchan
    FPEchx=squeeze(FPE_all(ik,:,:)); 
    SBCchx=squeeze(SBC_all(ik,:,:));
    FPE_mean_all(ik,:)=mean(FPEchx,1);
    SBC_mean_all(ik,:)=mean(SBCchx,1);
    FPE_mean_all=mean(FPE_all,2); 
    FPE_mean_all=squeeze(FPE_mean_all);
%     % Figure commented 8.7.2015
%     figure; plot(FPE_mean_all(ik,:));title(['mean FPE (blue) & SBC (red) of channel' num2str(ik)]); axis tight; hold on;
%     plot(SBC_mean_all(ik,:), 'r');xlabel('model order');
%     stemp=['FPE_SBC_mean_ch' num2str(ik)];
%     saveas(gcf, stemp, 'fig');
end
 clear C FPE FPE1 FPE1mean FPEchx N SBC SBCchx b2 ans ik jk k m n 
 clear num_epochs pend pstart 
 save SBC_FPE SBC_mean_all FPE_mean_all SBC_all FPE_all  % saves the results in the variable SBC_FPE 
 
 %% new add 15-11-2011
  for kk=1:nchan
     a1=min(FPE_mean_all(kk,:));
     p1=find(FPE_mean_all(kk,:)==a1);
     p1_all(kk)=p1;
  end
  
 % Figure commented 8.7.2015
  %figure; plot(p1_all,'*')

for kk=1:nchan
     a2=min(SBC_mean_all(kk,:));
     p2=find(SBC_mean_all(kk,:)==a2);
     p2_all(kk)=p2;
end

final_mean_p=mean(p2_all);
p=ceil(final_mean_p);
%  % Figure commented 8.7.2015
% figure; plot(p2_all,'*'); legend;
