close all
clear all
clc

%Import data 
filename = 'Data_sun_earth_mars.xlsx';
sheet = 4;
xlRange = 'A2:I7303';

[v,T,vT] = xlsread(filename, sheet, xlRange);
Sun_x=v(:,1); 
Sun_y=v(:,2); 
Sun_z=v(:,3); 
Earth_x=v(:,4); 
Earth_y=v(:,5); 
Earth_z=v(:,6); 
Mars_x=v(:,7); 
Mars_y=v(:,8); 
Mars_z=v(:,9); 


figure
scatter3(Sun_x,Sun_y,Sun_z,'k')
hold on
scatter3(Earth_x,Earth_y,Earth_z,'bx')
hold on 
scatter3(Mars_x,Mars_y,Mars_z,'rx')
legend('Sun','Earth','Mars','Fontsize',14)
xlabel('x','Fontsize',14)
ylabel('y','Fontsize',14)
zlabel('z','Fontsize',14)



