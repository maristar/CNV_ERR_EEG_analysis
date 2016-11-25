function [threshold] = plot2deeg3_frontal(connectivity_measure, XYZ, B)
% This program plots the location of EEG electrodes,inside a sphere,looking like a head fantom.
% The coordinates should be in cartesian coordinates and in a single variable
% xyz -- coordinate format X Y Z -we use only X and Y here
% the connectivity_measure must be in format N x N x number of repetitions
% where N is the number of channels
% B is cell with the names of the electrodes (so B has N elements as the number of channels).
% works Well 5-5-2011
% Oslo, 16-5-2011, Maria L. Stavrinou
nchan=length(XYZ);
% You can define the center of the sphere and the radius R. 
center2=[0 0];
% Radius can be the max of the radius from the coordinates 
% for k=1:length(xyz)
%     Radiusxyz(k)=(sqrt((xyz(k,1)^2)+(xyz(k,2)^2)+(xyz(k,3)^2)));
% end
% Radius=max(Radiusxyz);

% Or you can define the radius 
xyz=XYZ;
Radius=max(sqrt(XYZ(:,1).^2+XYZ(:,2).^2+XYZ(:,3).^2));
h=plot(-xyz(:,2), xyz(:,1), '+');
hold on; h=circle2(center2,Radius,100,'b-');

%make the nose
xtria=[-0.1*Radius 0.1*Radius 0]; ytria=[Radius Radius (Radius+(Radius/10))];
h=fill(xtria, ytria, 'w');
axis tight
% Put labels of the electrodes for that you need a cell variable containing
% the name of the electrodes in the correct order.
for k=1:length(B)
text(-xyz(k,2), xyz(k,1), B(k), 'FontSize',14, 'HorizontalAlignment', 'center','VerticalAlignment', 'bottom'); end

% make the arrow tip is the function ... draw_arrow_tip
% Plot MEAN connectivity by any measure

threshold=(3)*max(squeeze(max(squeeze(mean(connectivity_measure)))));


hold on;axis(axis)
for k=1:nchan
    for jj=1:nchan
        if jj~=k
            result_temp=squeeze(connectivity_measure(k,jj,:));
            mean_result_temp=mean(result_temp);
            if (mean_result_temp>(threshold/2)) && (mean_result_temp < threshold ) ;
                integerL=1;
                axis tight; axis(axis)
                line([-xyz(k,2) -xyz(jj,2)], [xyz(k,1) xyz(jj,1)], 'Marker', '.', 'LineWidth', integerL*mean_result_temp);
                draw_arrow_tip([-xyz(k,2) xyz(k,1)], [-xyz(jj,2) xyz(jj,1)], integerL/10); 
            elseif mean_result_temp > threshold, integerL=3;%0.16 a good value here
                axis tight; axis(axis)
                line([-xyz(k,2) -xyz(jj,2)], [xyz(k,1) xyz(jj,1)], 'Marker', '.', 'LineWidth', 3*integerL*mean_result_temp); % 'Color', [0 1 1]
                draw_arrow_tip([-xyz(k,2) xyz(k,1)], [-xyz(jj,2) xyz(jj,1)],integerL/10);%0.16 a good value here
            end
            hold on;
            clear result_temp
        end
    end
end


