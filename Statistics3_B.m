                                                                                                                                                                                                                      % Simple statistics
%% fyllo1 -coh, 
% list of names

% filenamesFolders={'kol27may08'};%PART 1 'kol27may08', ...
  %,'Sch21sep09', 'SRA09jun08', 'ped09sep08', 'hau14feb07','hol03mar08','jur08apr08', 'kol27may08', 'ped09sep08'
  %'ped09sep08', 'rau14aug07', 'Gun07sep09', 'hau14feb07','hol03mar08-raw', 'jur08apr08', 'new kle', 'Sch21sep09', 'SRA09jun08', 
  %STATSFINAL.(textmeasure).measure2=zeros(length(filenamesFolders), 2,4);
for jjk=1:length(filenamesFolders)
   disp(jjk)
    cd D:\DATA_MARIA
    cd(filenamesFolders{jjk})
    
filenamesFolders={'kol27may08'};%PAR

    % kane pragmata
load resultscor

% load the results from assymetry index giving the top 15 couples at sleep 
board=resultscor.(textmeasure).couple_conn;
board_values=resultscor.(textmeasure).couple_conn_values;
name=resultscor.filename(1:end-4);
resultscor.(textmeasure).chan_strength_norm=chan_strength_norm;
resultscor.(textmeasure).chan_appearance_norm=chan_appearance_norm;
% start the function --count appearance of single channel in the couples



%% GO

% START STATS
%% Measure 1: average connectivity of those 15 couples for day and night / per
% subject --> mannwhitney or ttest between day and night (repetitions are the subjects)
measure1_go=resultscor.(textmeasure).couple_conn_values;
measure1_nogo=resultscor.(textmeasure).couple_conn_values_nogo;
STATS.(textmeasure).measure1.description='average connectivity of those 15 couples for Go and NOGO --> mannwhitney or ttest between day and night (repetitions are the subjects)';
STATS.(textmeasure).measure1.go=mean(measure1_go);
STATS.(textmeasure).measure1.nogo=mean(measure1_nogo);
% how to add from other subjects 
%STATS.measure1.awake(2)=100;
allmeasure1=zeros(1,2);  % allmeasure1 nsubj x 2 (meansleep vs meanawake)
allmeasure1(1,1)=STATS.(textmeasure).measure1.go; 
allmeasure1(1,2)=STATS.(textmeasure).measure1.nogo;

clear measure1_awake measure1_sleep
STATS.(textmeasure).measure1=allmeasure1;
%% Measure 2: Spiking channels strength and appearance comparison between night and day.
STATS.(textmeasure).measure2.description='Test whether the strength of the spiking channel, changed significantly between night and day';
name=resultscor.filename(1:3);
switch lower(name)
case 'tba', spike_chan_index=5; % P3
case 'sra', spike_chan_index=[2,4]; % F4, C4
case 'kol', spike_chan_index=[5,7]; % P3, O1
case 'ped', spike_chan_index=[3,5]; % C3, P3 
case 'sch', spike_chan_index=5; % P3
case 'hau', spike_chan_index=[3,4];% C3 C4
    case 'hol', spike_chan_index=[7]; %Pz, O1
    case 'kle', spike_chan_index=[3,5]; % C3, P3 % added 1-8-2011
    case 'jur', spike_chan_index=[4,6];%C4 P4
    case 'rau', spike_chan_index=[9,10]; %P7 P8 
end

NumSpikes=length(spike_chan_index);

STATS.(textmeasure).measure2.spike_strength_sleep=resultscor.(textmeasure).chan_strength_norm_sleep(spike_chan_index);
STATS.(textmeasure).measure2.spike_strength_awake=resultscor.(textmeasure).chan_strength_norm_awake(spike_chan_index);

STATS.(textmeasure).measure2.spike_chan_appearance_sleep=resultscor.(textmeasure).chan_appearance_sleep(spike_chan_index);
STATS.(textmeasure).measure2.spike_chan_appearance_awake=resultscor.(textmeasure).chan_appearance_awake(spike_chan_index);

allmeasure2=zeros(NumSpikes,4);
allmeasure2(:,1)=resultscor.(textmeasure).chan_strength_norm_sleep(spike_chan_index);
allmeasure2(:,2)=resultscor.(textmeasure).chan_strength_norm_awake(spike_chan_index);
allmeasure2(:,3)=resultscor.(textmeasure).chan_appearance_sleep(spike_chan_index);
allmeasure2(:,4)=resultscor.(textmeasure).chan_appearance_awake(spike_chan_index);
STATS.(textmeasure).measure2.allmeasure2=allmeasure2;

%% Measure 3: 8 first couples during day, see their activity during night
STATS.(textmeasure).measure3.description='8 first couples during day, see their activity during night';

NumC=8;% Number of couples to do the calculation
%awake data: the 15 values of the 15 most strong couples and the 15 couples names, for example 'C4-P4'
a=resultscor.(textmeasure).couple_awake_values;
atext=resultscor.(textmeasure).couple_awake;
% sleep data: the same for sleep as above 
s=resultscor.(textmeasure).couple_sleep_values;
stext=resultscor.(textmeasure).couple_sleep;

%% WE need the 8 first couples during the day and then to see the same
% couples during the night.
matrix1_day=zeros(NumC, 2);
matrix1_day_text=cell(NumC, 1);
for jj=1:NumC;
    matrix1_day_text{jj}=atext{jj};
    matrix1_day(jj,1)=a(jj);
    k1=strmatch(atext{jj}, stext);
    if k1~=0 
        matrix1_day(jj,2)=s(k1);
    else
        matrix1_day(jj,2)=0;
    end
clear k1
end
STATS.(textmeasure).measure3.matrix_day=matrix1_day;
STATS.(textmeasure).measure3.matrix_day_text=matrix1_day_text;


%% Here we need the 8 first couples during the night and then to see the
% same duting the day 
matrix2_night=zeros(NumC, 2);
matrix2_night_text=cell(NumC,1);
for jj=1:NumC
    matrix2_night(jj,1)=s(jj); % 8 x 2 double
    matrix2_night_text{jj}=stext{jj};  % 1x 8 cell
    k1=strmatch(stext{jj}, atext);
    if k1~=0 
        matrix2_night(jj,2)=a(k1);
    else
        matrix2_night(jj,2)=0;
    end
    clear k1
end

STATS.(textmeasure).measure3.matrix2_night=matrix2_night;
STATS.(textmeasure).measure3.matrix2_night_text=matrix2_night_text;

% create allmeasure3
allmeasure3=[matrix1_day matrix2_night];
STATS.(textmeasure).measure3.allmeasure3=allmeasure3;

clear matrix2_night matrix2_night_text jj k1 atext a s NumC
clear matrix1_day matrix1_day_text jj k1 atext stext
% cd(resultscor.pathname)
% cd ..
% cd RAW

% End Measure 3
%% MEASURE 4
% compare the percentages of the spiking channels with the percentage of
% the most strong and the most appearances 
% We need, parts of the measure 2. 

%% STRENGTH
% Strength of the spiking channels and strength of the more strong channel.
% strength of spiking channel at sleep and awake
STATS.(textmeasure).measure4.spike_strength_sleep=resultscor.(textmeasure).chan_strength_norm_sleep(spike_chan_index);
STATS.(textmeasure).measure4.spike_strength_awake=resultscor.(textmeasure).chan_strength_norm_awake(spike_chan_index);

%sleep strength max
temp1=resultscor.(textmeasure).chan_strength_norm_sleep;
[maxchan, ind]=max(temp1);
STATS.(textmeasure).measure4.strength_sleep_max=temp1(ind);
clear ind
clear maxchan
clear temp1

% awake strength max
temp1=resultscor.(textmeasure).chan_strength_norm_awake;
[maxchan, ind]=max(temp1);
STATS.(textmeasure).measure4.strength_awake_max=temp1(ind);
clear ind
clear maxchan
clear temp1

%% APPEARANCE RATE
% percentage of appearance max 
STATS.(textmeasure).measure4.spike_chan_appearance_sleep=resultscor.(textmeasure).chan_appearance_sleep(spike_chan_index);
STATS.(textmeasure).measure4.spike_chan_appearance_awake=resultscor.(textmeasure).chan_appearance_awake(spike_chan_index);

%sleep appearance max 
temp1=resultscor.(textmeasure).chan_appearance_sleep;
[maxchan, ind]=max(temp1);
STATS.(textmeasure).measure4.appearance_sleep_max=temp1(ind);
clear ind
clear maxchan
clear temp1


% awake appearance max
temp1=resultscor.(textmeasure).chan_appearance_awake;
[maxchan, ind]=max(temp1);
STATS.(textmeasure).measure4.appearance_awake_max=temp1(ind);
clear ind
clear maxchan
clear temp1

% a variable with all inside named allmeasure4
allmeasure4=zeros(NumSpikes,8);%allmeasure4=zeros(22,8)
% strength 
%1)spike sleep 2) max channel sleep 3) spike awake 4) max channel awake
allmeasure4(:,1)=STATS.(textmeasure).measure4.spike_strength_sleep;
allmeasure4(:,2)=STATS.(textmeasure).measure4.strength_sleep_max;
allmeasure4(:,3)=STATS.(textmeasure).measure4.spike_strength_awake;
allmeasure4(:,4)=STATS.(textmeasure).measure4.strength_awake_max;

% appearance 
%5)spike sleep 6) max channel sleep 7) spike awake 8) max appear channel awake
allmeasure4(:,5)=STATS.(textmeasure).measure4.spike_chan_appearance_sleep;
allmeasure4(:,6)=STATS.(textmeasure).measure4.appearance_sleep_max;
allmeasure4(:,7)=STATS.(textmeasure).measure4.spike_chan_appearance_awake;
allmeasure4(:,8)=STATS.(textmeasure).measure4.appearance_awake_max;
STATS.(textmeasure).measure4.allmeasure4=allmeasure4;
save STATS STATS
% [p,h]=ranksum(allmeasure4(:,1), allmeasure4(:,2)) ; % statistical
% significnace between spiking electrode strength and max channel strength
% at Sleep
%[p,h]=ranksum(allmeasure4(:,3), allmeasure4(:,4)) ; % statistical
% significnace between spiking electrode strength and max channel strength
% at awake
%[p,h]=ranksum(allmeasure4(:,5), allmeasure4(:,6)) ; % statistical
% significnace between spiking electrode appearance rate and max channel strength
% at sleep
%[p,h]=ranksum(allmeasure4(:,7), allmeasure4(:,8)) ; % statistical
% significnace between spiking electrode appearance rate and max channel strength
% at awake
%STATS.measure4.pvalues=
% 
%% MEASURE 5
% Strength and appearance rate of occipital electrodes. The occipital
% electrodes are 7 and 8 in index files for raw. They are in 7, 8 in
% Laplacian 
% O1 and o2 strength go together. and compare 01 and 02 during day and
% night

% O1 strength sleep and awake
STATS.(textmeasure).measure5.O1_strength_sleep=resultscor.(textmeasure).chan_strength_norm_sleep(7);
STATS.(textmeasure).measure5.O1_strength_awake=resultscor.(textmeasure).chan_strength_norm_awake(7);

%O1 appearance 
STATS.(textmeasure).measure5.O1_appearance_sleep=resultscor.(textmeasure).chan_appearance_sleep(7);
STATS.(textmeasure).measure5.O1_appearance_awake=resultscor.(textmeasure).chan_appearance_awake(7);

% O2 strength sleep and awake
STATS.(textmeasure).measure5.O2_strength_sleep=resultscor.(textmeasure).chan_strength_norm_sleep(8);
STATS.(textmeasure).measure5.O2_strength_awake=resultscor.(textmeasure).chan_strength_norm_awake(8);

%O2 appearance 
STATS.(textmeasure).measure5.O2_appearance_sleep=resultscor.(textmeasure).chan_appearance_sleep(8);
STATS.(textmeasure).measure5.O2_appearance_awake=resultscor.(textmeasure).chan_appearance_awake(8);


% %1)occipital strength sleep 2) occipital strength awake 
% allmeasure5_sleep_occipital(:,1)=resultscor.(textmeasure).chan_strength_norm_sleep(7);
% allmeasure5_sleep_occipital(:,2)=resultscor.(textmeasure).chan_strength_norm_sleep(8);
% allmeasure5_awake_occipital(:,1)=resultscor.(textmeasure).chan_strength_norm_awake(7);
% allmeasure5_awake_occipital(:,2)=resultscor.(textmeasure).chan_strength_norm_awake(8);
% 
% % 1) occipital appearance sleep 2) occipital appearance awake
% allmeasure5_sleep_occipital_ap(:,1)=resultscor.(textmeasure).chan_appearance_sleep(7);
% allmeasure5_sleep_occipital_ap(:,2)=resultscor.(textmeasure).chan_appearance_sleep(8);
% allmeasure5_awake_occipital_ap(:,1)=resultscor.(textmeasure).chan_appearance_awake(7);
% allmeasure5_awake_occipital_ap(:,2)=resultscor.(textmeasure).chan_appearance_awake(8);
% allmeasure5
allmeasure5=zeros(2,4); % allmeasure5=zeros(NumData*2,4)
allmeasure5(1,1)=resultscor.(textmeasure).chan_strength_norm_sleep(7);
allmeasure5(2,1)=resultscor.(textmeasure).chan_strength_norm_sleep(8);
allmeasure5(1,2)=resultscor.(textmeasure).chan_strength_norm_awake(7);
allmeasure5(2,2)=resultscor.(textmeasure).chan_strength_norm_awake(8);

% 1) occipital appearance sleep 2) occipital appearance awake
allmeasure5(1,3)=resultscor.(textmeasure).chan_appearance_sleep(7);
allmeasure5(1,4)=resultscor.(textmeasure).chan_appearance_sleep(8);
allmeasure5(2,3)=resultscor.(textmeasure).chan_appearance_awake(7);
allmeasure5(2,4)=resultscor.(textmeasure).chan_appearance_awake(8);
%
STATS.(textmeasure).measure5.allmeasure5=allmeasure5;
save STATS STATS
%end of measure 5
%% 

%% Save the results into resultscor
save STATS STATS
resultscor.(textmeasure).stats=STATS;
save resultscor resultscor -v7.3
clear STATS

cd('D:\DATA_MARIA\')
STATSFINAL_lapl(jjk).(textmeasure).measure1(:,:)=allmeasure1;
STATSFINAL_lapl(jjk).(textmeasure).measure2(:,:)=allmeasure2;
% if NumSpikes==1
%     STATSFINAL2(jjk).(textmeasure).measure2(1,:)=allmeasure2;
%     STATSFINAL2(jjk).(textmeasure).measure2(2,:)=0;
% end
%STATSFINAL2(jjk).(textmeasure).measure2(:,:)=allmeasure2;
STATSFINAL_lapl(jjk).(textmeasure).measure3(:,:)=allmeasure3;
STATSFINAL_lapl(jjk).(textmeasure).measure4(:,:,:,:,:,:,:,:)=allmeasure4;
STATSFINAL_lapl(jjk).(textmeasure).measure5(:,:,:,:)=allmeasure5;
% STATSFINAL_LAPL.(textmeasure).measure1(jjk,:,:)=allmeasure1;
% % if NumSpikes==2
% %      STATSFINAL_LAPL.(textmeasure).measure2(jjk, 1,:)=allmeasure2;
% %      STATSFINAL_LAPL.(textmeasure).measure2(jjk, 2,:)=0;
% %  end
% STATSFINAL.(textmeasure).measure2(jjk, :,:)=allmeasure2;
% STATSFINAL_LAPL.(textmeasure).measure3(jjk,:,:)=allmeasure3;
% STATSFINAL_LAPL.(textmeasure).measure4(jjk,:,:,:,:,:,:,:,:)=allmeasure4;
% STATSFINAL_LAPL.(textmeasure).measure5(jjk,:,:,:,:)=allmeasure5;
% 
save STATSFINAL STATSFINAL
% clear all
clear NumSpikes STATS resultscor allmeasure1 allmeasure2 allmeasure3 allmeasure4 allmeasure5 spike_chan_index
    %MEASURE1_AWAKE(jjk,:)=resultscor.(textmeasure).stats.measure1;
end

%filenamesFolders={'hau14feb07','hol03mar08','jur08apr08', 'kol27may08', 'ped09sep08','rau14aug07', 'Sch21sep09','SRA09jun08','tba10jan08'};
filenamesFolders={'kle03aug09b'};
