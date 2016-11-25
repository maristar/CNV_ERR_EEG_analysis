% select the appropriate model order 
NN=size(data_1, 1)
for kk=1:nchan
for nn=1:55
   data_1=data(kk,:,nn)';
    for p=1:20
    [popt w, A, C, sbc, fpe, th]=arfit_m(data_1, p, p);
    BIC=log(C)+p/(NN*log(NN));
    BIC_all(kk,p,nn,:)=BIC;
    fpe_all(kk,p,nn,:)=fpe;
    end
end
end

% how we see the fpe
fpe=squeeze(fpe_all(1,:,:));
for kk=1:26; figure(3); plot(fpe(:,kk)); hold on; end



for p=1:40 
        [ar_coeffs,NoiseVariance, reflect_coeffs]=aryule(interpol_intervaloRR, p);
        %figure(100); plot(reflect_coeffs); hold on; 
        %reflect_coeffs_all(q,p,;)=reflect_coeffs;
        Noise_variance_all(q,p,:)=NoiseVariance;
        BIC=log(NoiseVariance)+p/(NN*log(NN));
        BIC_all(q,p,:)=BIC
        BIC=[];
        NoiseVariance=[];
        reflect_coeffs=[];
       % figure; plot(mean(BIC_all(:,1)))
end