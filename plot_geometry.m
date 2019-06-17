x=load('../matlabx.out');
y=load('../matlaby.out');
z=load('../matlabz.out');


hold off;


h=surf(x,y,z);
set(h,'FaceLighting','gouraud','FaceColor',[0.65 0.65 0.65],'EdgeColor','none',...
        'SpecularColorReflectance',0.4,...
        'SpecularStrength',0.2,'DiffuseStrength',0.9,'AmbientStrength',0.8);
light('Position',[1 -2 1]);
axis equal;
material dull;
view(-10,30);
set(get(gca,'XLabel'),'String','x');
set(get(gca,'YLabel'),'String','y');
set(get(gca,'ZLabel'),'String','z');
hold on;

% Coordinate Axes: 
% To see the standard orientation of the particle (all Euler angles
% equal to zero), it is best to plot the particle WITH
% the coordinate axes. For presentation purposes, it is best 
% to plot the particle WITHOUT the coordinate
% system.
% To remove the coordinate axes, uncomment the next line:

%axis vis3d off;
