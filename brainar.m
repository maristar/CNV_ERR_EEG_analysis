function [br_areas]=brainar(chan_strength_norm, chan_appearance_norm)
% we use it after the strength_appearance function in order to calculate
% the average connectivity of each brain area.

CS=chan_strength_norm;
CA=chan_appearance_norm;
% Frontal Left
FL=(CS(7)+CS(9)+CS(8))/3;

% Frontal right
FR=(CS(1)+CS(28)+CS(29))/3;

% Frontal Z
FZ=(CS(2)+CS(5)+CS(6))/3;

% Central Z
CZ=(CS(3)+CS(25)+CS(4))/3;

%Central Right
CR=(CS(24)+CS(27)+CS(26))/3;

%Central Left
CL=(CS(10)+CS(11)+CS(12))/3;
% Parietal Z
PZ=(CS(16)+CS(18))/2;

% Parietal Left
PL=(CS(13)+CS(14)+CS(15))/3;

% Parietal Right
PR=(CS(21)+CS(22)+CS(23))/3;

% Occipital
OZ=(CS(17)+CS(19)+CS(20))/3;
br_areas.FL=FL;
br_areas.FR=FR;
br_areas.CZ=CZ;
br_areas.CL=CL;
br_areas.CR=CR;
br_areas.PZ=PZ;
br_areas.PL=PL;
br_areas.PR=PR;
br_areas.OZ=OZ;
br_areas.FZ=FZ;


