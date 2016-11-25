function [grandaverage, XYZ]=grandaverager_m(DOI, trigger, textmeasure, nchan)
grandaverage=zeros(nchan);
cd(DOI)
cd(upper(trigger))
    files=dir('*');
    for kk=1:length(files); 
        filenames{kk,:}=files(kk,:).name;
    end
    ND=kk; clear kk;
    meantemp=0;
for kkm=3:ND
        cd(DOI)
        cd(upper(trigger))%% sos to make upcase
        %disp(filenames(kkm));
        disp(kkm-2)
        stemp=([filenames{kkm}]);
        cd(stemp)
        load resultscor
        XYZ=resultscor.XYZ;
        %resultsALL.s=resultscor.s;
        switch textmeasure
            case 'pcor', meantemp=squeeze(mean(resultscor.(textmeasure).result2,3));
            case 'cor', meantemp=squeeze(mean(resultscor.(textmeasure).result1,3));
            case {'DTFdelta','DTFtheta','DTFalpha','DTFbeta','DTFgamma'}, meantemp=resultscor.resultsDTF.(textmeasure).meangamma; %for %dtf
        end
        grandaverage=meantemp+grandaverage;
        clear meantemp resultscor
end
    grandaverage=grandaverage./(ND-2);
    
    XYZ=XYZ;