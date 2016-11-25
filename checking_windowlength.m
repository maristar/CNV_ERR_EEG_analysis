% testing window length
win=100 % ms
windp=floor(win*Fs/1000)
   
%% try 2
    NN=size(st1,2);
    windp=floor(win*Fs/1000);
    num_ep=floor(NN/win);
    winwind=hamming(windp+1);
for k=1:num_ep;
    temp1=st1((k+(k*windp)-windp):(k+((k+1)*windp)-windp))';
    temp2=st2((k+(k*windp)-windp):(k+((k+1)*windp)-windp))';
    wtemp1=convn(winwind, temp1);
    wtemp2=convn(winwind, temp2);
    x7(k,:,:)=corr(temp1,temp2);
end
xcorre=xcorr(st1, st2);
timeVec7=(1:(length(x7))).*length(xcorre)/length(x7);
% timeVectorSI=(1:length(vSI_norm)).*600/3600;


figure; subplot(2,1,1); plot(xcorre); subplot(2,1,2); plot(timeVec7, x7, 'r*'); 

 

