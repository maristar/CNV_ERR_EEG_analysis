   
cd(DOI)
       
    s=N;
    h=0;
    clear index_pool kj kl 
    index_pool={1:nchan*nchan}; %% to change this to the value of possible combinations
    for kl=1:nchan
        for kj=1:nchan
            if (kl~=kj)
                h=h+1;
                index_pool{h,:}=[num2str(kj) num2str(kl)];
                a=strcmp(index_pool, [num2str(kl),num2str(kj)]);
                onesz=find(a);
                if isempty(onesz)
                    go=squeeze(grandaverage(kj,kl,:));   
                    nogo=squeeze(grandaverage_nogo(kj,kl,:));
                    %[p,h] = signrank(go, nogo);
                    [h1,p,ci] =ttest2(go,nogo,0.05,'both','unequal');
                    significance(kl,kj,:)=p;
                    h_value(kl,kj,:)=h1;
                    clear p, h1
                end % if
                clear onesz a 
            end % if (kl~=kj)
        end %kj
    end %kl
    for jj=1:nchan;
        h_value(jj,jj)=0;
        significance(jj,jj)=0;
    end
 
    % a way to plot and see what is significant and what is not. 
    for kl=1:nchan
      for kj=1:nchan
         if significance(kl,kj)<0.05
             significance_plot(kl,kj)=1;
            else
             significance_plot(kl,kj)=0;
         end
      end
    end
    resultsALL.(textmeasure).ttest=significance;
    resultsALL.(textmeasure).ttest_plot=significance_plot;
    resultsALL.(textmeasure).ttest_h=h_value;
    cd(textmeasure)
    save resultsALL resultsALL -v7.3
       
    figure; imagesc(h_value);title([textmeasure 'h value']);
    set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan); 
    set(gca, 'YTickLabel', s); set(gca, 'XTickLabel', s);
    axis xy; axis tight; colorbar('location','EastOutside')
    
    figure; imagesc(significance_plot);title([textmeasure 'p value']);
    set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan); 
    set(gca, 'YTickLabel', s); set(gca, 'XTickLabel', s);
    axis xy; axis tight; colorbar('location','EastOutside')
    title('significance p value')
    
    clear grandaverage_nogo grandaverage
end