function plot2dhead(connectivity_measure, xyz, B)
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

% Or you can define the radius 

Radius=0.0999;
figure; h=plot(xyz(:,2), xyz(:,1), '+');
hold on; h=circle2(center2,Radius,100,'b-');

%make the nose
xtria=[-0.01 0.01 0]; ytria=[0.1 0.1 0.11];
h=fill(xtria, ytria, 'w');
axis tight
% Put labels of the electrodes for that you need a cell variable containing
% the name of the electrodes in the correct order.
for k=1:length(B)
text(xyz(k,2), xyz(k,1), B(k), 'FontSize',14, 'HorizontalAlignment', 'center','VerticalAlignment', 'bottom'); end


% Plot MEAN connectivity by any measure
hold on;
for k=1:10
    for jj=1:10
        if jj~=k
            result_temp=squeeze(connectivity_measure(k,jj,:));
            mean_result_temp=mean(result_temp);
            if mean_result_temp>median(connectivity_measure)/2, integerL=3;
            else integerL=1;
            end
            %if mean_result_temp<0.2; integerL=0.01; end
            %arrow([xyz(k,2) xyz(jj,2)], [xyz(k,1) xyz(jj,1)], 'Marker', '.', 'LineWidth', integerL)
            h=line([xyz(k,2) xyz(jj,2)], [xyz(k,1) xyz(jj,1)], 'Marker', '.', 'LineWidth', integerL);
            %vectarrow([xyz(k,2) xyz(jj,2)], [xyz(k,1) xyz(jj,1)])
           % draw_arrow([xyz(k,2) xyz(jj,2)], [xyz(k,1) xyz(jj,1)], 0.08)
%             h=arrow([xyz(k,2) xyz(jj,2)], [xyz(k,1) xyz(jj,1)],26,'BaseAngle',60);
%             set(h, 'Marker', '.', 'LineWidth', integerL);
            hold on;
            clear result_temp
        end
    end
end


