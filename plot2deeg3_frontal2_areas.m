function plot2deeg3_frontal2_areas(connectivity_measure, XYZ, numLeft, Far, Car, Par, thr)
% This program plots the location of EEG electrodes,inside a sphere,looking like a head fantom.
% The coordinates should be in cartesian coordinates and in a single variable
% xyz -- coordinate format X Y Z -we use only X and Y here
% the connectivity_measure must be in format N x N x number of repetitions
% where N is the number of channels
% B is cell with the names of the electrodes (so B has N elements as the number of channels).
% works Well 5-5-2011
% Oslo, 16-5-2011, Maria L. Stavrinou
% for use with intrahemispheric -16-12-2011 MLS num is numRight or numLeft

% You can define the center of the sphere and the radius R. 
center2=[0 0];

XYZareaL=XYZ(numLeft,:);
xyz=XYZareaL;

Radius=max(sqrt(xyz(:,1).^2+xyz(:,2).^2+xyz(:,3).^2));
% h=plot(-xyz(:,2), xyz(:,1), '+');
% hold on; 

h=circle2(center2,Radius,100,'b-');
hold on;

FL1=(xyz(Far(1),2)+xyz(Far(2),2)+xyz(Far(3),2))/3; 
disp('ok fl1')
FL2=(xyz(Far(1),1)+xyz(Far(2),1)+xyz(Far(3),1))/3;
plot(-FL1, FL2, 'r*');
FLx=[-FL1 FL2];

CL1=(xyz(Car(1),2)+xyz(Car(2),2)+xyz(Car(3),2))/3; CL2=(xyz(Car(1),1)+xyz(Car(2),1)+xyz(Car(3),1))/3;
plot(-CL1, CL2, 'r*');
CLx=[-CL1 CL2];

PL1=(xyz(Par(1),2)+xyz(Par(2),2)+xyz(Par(3),2))/3; PL2=(xyz(Par(1),1)+xyz(Par(2),1)+xyz(Par(3),1))/3;
plot(-PL1, PL2, 'r*');
PLx=[-PL1 PL2];

AR=[FLx; CLx; PLx];

%B={'FL' 'CL' 'PL'}
%make the nose
xtria=[-0.1*Radius 0.1*Radius 0]; ytria=[Radius Radius (Radius+(Radius/10))];
h=fill(xtria, ytria, 'w');
axis tight
% Put labels of the electrodes for that you need a cell variable containing
% the name of the electrodes in the correct order.
%  for k=1:length(AR)
%  text(-AR(k,2), AR(k,1), B(k), 'FontSize',10, 'HorizontalAlignment', 'center','VerticalAlignment', 'bottom'); end

% make the arrow tip is the function ... draw_arrow_tip
% Plot MEAN connectivity by any measure

threshold=thr; %(3)*max(squeeze(max(squeeze(mean(connectivity_measure(num))))));

nchan=length(AR);
hold on;axis(axis)
for k=1:nchan
    for jj=1:nchan
        if jj~=k
            result_temp=squeeze(connectivity_measure(k,jj,:));
            mean_result_temp=mean(result_temp);
            if (mean_result_temp>(threshold/2)) && (mean_result_temp < threshold ) ;
                integerL=1;
                axis tight; axis(axis)
                line([AR(k,1) AR(jj,1)], [AR(k,2) AR(jj,2)], 'Marker', '.', 'LineWidth', integerL);
                draw_arrow_tip([AR(k,1) AR(k,2)], [AR(jj,1) AR(jj,2)], integerL/10); 
            elseif mean_result_temp > threshold, 
                integerL=3;
%                 disp('mean res')
%                 disp(mean_result_temp)
%                 disp('integerL')
%                 disp(30*integerL*mean_result_temp)%0.16 a good value here
                axis tight; axis(axis)
                line([AR(k,1) AR(jj,1)], [AR(k,2) AR(jj,2)], 'Marker', '.', 'LineWidth', integerL);
                draw_arrow_tip([AR(k,1) AR(k,2)], [AR(jj,1) AR(jj,2)], integerL/10); 
            end
            hold on;
            clear result_temp
        end
    end
end


