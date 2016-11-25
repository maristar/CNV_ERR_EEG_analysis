function plot2dhead_frontal2(connectivity_measure, XYZ, B, thr)
% 12-12-2011 for frontal lobe analysis for grandaverages
% THIS IS THE CORRECT FUNCTION
% This program plots the location of EEG electrodes,inside a sphere,looking like a head fantom.
% The coordinates should be in cartesian format and in a single variable
% xyz -- coordinate format X Y Z -we use only X and Y here
% the connectivity_measure must be in format N x N x number of repetitions 
% B is cell with the names of the electrodes (so B has N elements as the number of electrodes).
% works Well 5-5-2011
% You can define the center of the sphere and the radius R. 
center2=[0 0];
% Radius can be the max of the radius from the coordinates 
% for k=1:length(xyz)
%     Radiusxyz(k)=(sqrt((xyz(k,1)^2)+(xyz(k,2)^2)+(xyz(k,3)^2)));
% end
% Radius=max(Radiusxyz);
nchan=length(B);
% Or you can define the radius 

Radius=max(sqrt(XYZ(:,1).^2+XYZ(:,2).^2+XYZ(:,3).^2));
%figure; 
h=plot(-XYZ(:,2), XYZ(:,1), '+');
% figure; h=plot3(-XYZ(:,2), XYZ(:,1), XYZ(:,3), '+');for k=1:length(B)
% text(-XYZ(k,2), XYZ(k,1), XYZ(k,3), B(k), 'FontSize',12, 'HorizontalAlignment', 'center','VerticalAlignment', 'bottom'); end
hold on; h=circle2(center2,Radius,100,'b-');

%make the nose
xtria=[-0.1*Radius 0.1*Radius 0]; ytria=[Radius Radius (Radius+(Radius/10))];
h=fill(xtria, ytria, 'w');
axis tight
% Put labels of the electrodes for that you need a cell variable containing
% the name of the electrodes in the correct order.
for k=1:length(B)
text(-XYZ(k,2), XYZ(k,1), B(k), 'FontSize',10, 'HorizontalAlignment', 'center','VerticalAlignment', 'bottom'); end

% threshold=(3/2)*(max(squeeze(mean(connectivity_measure))));
threshold=thr;
% Plot MEAN connectivity by any measure
hold on;
for k=1:nchan
    for jj=1:nchan
        if jj~=k
            result_temp=squeeze(connectivity_measure(k,jj));
            mean_result_temp=mean(result_temp);
            if mean_result_temp>threshold;
                integerL=3;% mean(connectivity_measure)
                h=line([-XYZ(k,2) -XYZ(jj,2)], [XYZ(k,1) XYZ(jj,1)], 'Marker', '.', 'LineWidth', integerL);
            else integerL=1;
                %h=line([-XYZ(k,2) -XYZ(jj,2)], [XYZ(k,1) XYZ(jj,1)], 'Marker', '.', 'LineWidth', integerL/20);
            end
            %if mean_result_temp<0.2; integerL=0.01; end
            %arrow([xyz(k,2) xyz(jj,2)], [xyz(k,1) xyz(jj,1)], 'Marker', '.', 'LineWidth', integerL)
            %h=line([XYZ(k,2) XYZ(jj,2)], [XYZ(k,1) XYZ(jj,1)], 'Marker', '.', 'LineWidth', integerL/20);
            %vectarrow([xyz(k,2) xyz(jj,2)], [xyz(k,1) xyz(jj,1)])
           % draw_arrow([xyz(k,2) xyz(jj,2)], [xyz(k,1) xyz(jj,1)], 0.08)
%             h=arrow([xyz(k,2) xyz(jj,2)], [xyz(k,1) xyz(jj,1)],26,'BaseAngle',60);
%             set(h, 'Marker', '.', 'LineWidth', integerL);
            hold on;
            clear result_temp
        end
    end
end


