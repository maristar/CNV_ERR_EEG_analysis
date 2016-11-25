% 15-6-2011 follows the assymetry index for calculation of strength and
% rate of appearance of electrodes from the 15 more strong couples
function [chan_strength_norm, chan_appearance_norm]=strength_appearance(board, board_values,s, crank)
N=length(board);
for hh=1:length(s) % hh einai to kathe kanali
    k=strfind(board, s{hh}); %tha dwsei kati san  [] [] [] [] [4] [] [] [] [] [] [] [] [] [] [4];
    countemp=0; 
    tempstrength=0;
    for jj=1:N, a=strmatch(k{jj}, []); 
        if a==1,countemp=countemp+1; % otan einai keno [], dinei 1 
            tempstrength(jj)=0;
        else
            tempstrength(jj)=board_values(jj);
        end, 
        countappear=crank-countemp; 
%         disp(countappear)
        chan_countappear(hh)=countappear;
        chan_strength(hh)=sum(tempstrength);
    end
    
    clear k a countappear countemp tempstrength
end
total_chan_strength=sum(chan_strength(:));
total_chan_appearance=sum(chan_countappear(:));
chan_appearance_norm=100*chan_countappear/total_chan_appearance;%% this was changed 27-6-2011
chan_strength_norm=100*chan_strength/total_chan_strength; 