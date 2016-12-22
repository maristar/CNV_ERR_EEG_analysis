function [CON]=majorlist(average_conn, textmeasure, now, name, nchan, s, crank)
% DEFINE LIST WITH 10 or 15 (Crank) MOST ACTIVE ELECTRODE COUPLES 

average_conn_l=tril(average_conn,-1); % LOWER LEFT diagonal of array- since  partial correlation is symmetrical
average_conn_u=triu(average_conn,0); 

%crank=15; % the top channels (how many we want)
Megisti=average_conn;

%[Megisti,ind]=sort(average_conn(:),'descend');
% [r, c, v]=find(Megisti);% [row, columns, logicalvalue]=find non zeros in Megisti. Correction 2012 for elimination of zeros
% Megisti=Megisti(r);
% 
% Megisti = Megisti(1:end); ind=ind(1:end);
ind=1:length(Megisti);
[row col] = ind2sub(size(average_conn),ind);
couple_conn={};

% create the list with the most active electrode couples
for k=1:length(Megisti), couple_conn{k}=[s{row(k)} '-' s{col(k)}]; end;
CON.couple_conn=couple_conn;

% this is the values of the connectivity or correlation coefficient for the most active
% channels
for kk=1:length(Megisti), couple_conn_values(kk)=Megisti(kk); end;
CON.couple_conn_values=couple_conn_values;

% Send the results to excel file 
stempp=['MostCouples-' textmeasure]; % do not delete this!!!!
xlswrite(stempp, {name}, 'Sheet1', 'A1:A1')
titles={'most strong couples'};
xlswrite(stempp, (titles), 'Sheet1', 'A2:A2')
xlswrite(stempp,couple_conn, 'Sheet1', 'B2')
xlswrite(stempp,couple_conn_values, 'Sheet1', 'B3')

