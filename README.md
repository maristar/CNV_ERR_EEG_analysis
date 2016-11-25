# Frontal_lobe
CNV Connectivity analysis

This repository contains the programs for the connectivity analysis of patient data with frontal lobe lessions. 
The preprocessing was done using functions from EEGLAB and custom made functions from EEGLAB.
The clean and preprocessed datasets were then put into the connectivity analysis. 

The connectivity analysis on the EEG datasets was done with the Directed transfer Function (DTF) method. Especially a 
modified version of the original DTF algorithm was made that calculates the connectivity in short time intervals (see Ding et al. 2001)
This modifications were made by Maria Stavrinou on functions of the e-connectome software. (http://econnectome.umn.edu/).

The analysis was first run on electrodes and the connectivity between them was calculated. 
Then the connectivity between preselected areas consisting of about 3 electrodes each was calculated.  

These areas are the following:
Frontal left (FL)
Frontal Right(FR)
Central Left (CL)
Central Right (CR)
Parietal Left (PL)
Parietal Right (PR)


The programs give out excel files with values of connectivity, that are to be used for the statistical analysis. 
