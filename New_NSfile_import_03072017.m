clear all;
close all;
% Raw_path='E:\CNV KONNEKTIVITET_2016 fortsettelse av prosjekt fra 2012\R?data Front 1 Kopiert April 2017\Kontroller\';
% Analyzed_path='E:\CNV KONNEKTIVITET_2016 fortsettelse av prosjekt fra 2012\R?data Front 1 Kopiert April 2017\Kontroller\Analyzed_data\';
% Raw_path='E:\CNV KONNEKTIVITET_2016 fortsettelse av prosjekt fra 2012\R?data Front 1 Kopiert April 2017\Laterale\';
% Analyzed_path='E:\CNV KONNEKTIVITET_2016 fortsettelse av prosjekt fra 2012\R?data Front 1 Kopiert April 2017\Laterale\Analyzed_data\';
Raw_path='/Users/mstavrin/Documents/MATLAB/CNV_July_2017/'
% 'E:\CNV KONNEKTIVITET_2016 fortsettelse av prosjekt fra 2012\R?data Front 1 Kopiert April 2017\orbital\';
Analyzed_path='/Users/mstavrin/Documents/MATLAB/CNV_July_2017/Imported_data/'
'/Users/mstavrin/Documents/MATLAB/CNV_July_2017/'
% 'E:\CNV KONNEKTIVITET_2016 fortsettelse av prosjekt fra 2012\R?data Front 1 Kopiert April 2017\orbital\Analyzed_data\';

% Make a directory where the data will be saved. 
cd(Analyzed_path)
save_data_folder='Imported_data';
mkdir(save_data_folder)
saving_path=[Analyzed_path save_data_folder];

% Find the raw files we want to open
cd(Raw_path)
listing_raw=dir('*CNV*.raw');
Num_folders=length(listing_raw);
for kk=1:Num_folders
    names{kk,:}=listing_raw(kk).name;
end
clear kk

for kk=1:Num_folders
    Folder_name=names{kk,:};
    Folder_name_char=char(Folder_name);
    Raw_path_folder=[Raw_path Folder_name_char];
        
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_readegi(Raw_path_folder, [],[],'auto'); %[Raw_path Folder_name]
    temp_setname1=[Folder_name_char(1:end-18) '_imp'];
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',temp_setname1,'gui','off'); 
    EEG = eeg_checkset( EEG );
    eeglab redraw
    
    % Average reference
    EEG = pop_reref( EEG, [],'refloc',struct('labels',{'Cz'},'Y',{0},'X',{5.2239e-16},'Z',{8.5313},'sph_theta',{0},'sph_phi',{90},'sph_radius',{8.5313},'theta',{0},'radius',{0},'type',{'REF'},'ref',{''},'urchan',{[]},'datachan',{0}));
    temp_setname2=[temp_setname1 '_reref'];
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'savenew',temp_setname2,'gui','off'); 
    EEG = eeg_checkset( EEG );
    eeglab redraw
    
    % Load the correct locations 
    EEG=pop_chanedit(EEG, 'load',[],'load',{'/Users/mstavrin/Documents/MATLAB/CNV_July_2017/EGI129-HydroCel-Original.spf' 'filetype' 'sfp'});
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'savenew',temp_setname2,'gui','off'); 
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    eeglab redraw
    
    % We save the final set 
    cd(saving_path)
    %temp_setname2=[Folder_name_char(1:end-18) '_imp_reref'];
    temp_savename=[temp_setname2 '.set'];
    EEG=pop_saveset(EEG, 'filename', temp_savename);
    EEG = eeg_checkset( EEG );
    eeglab redraw
end
 