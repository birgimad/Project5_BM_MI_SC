close all
clear all
clc

filename = 'Test_2body_system.xlsx';
sheet = 1;
xlRange = 'E4:F14';

[v,T,vT] = xlsread(filename, sheet, xlRange);
time_step_RK4_earth=v(:,1); 
distance_RK4_earth=v(:,2);

filename = 'Test_2body_system.xlsx';
sheet = 1;
xlRange = 'I4:J14';

[v,T,vT] = xlsread(filename, sheet, xlRange);
time_step_VV_earth=v(:,1); 
distance_VV_earth=v(:,2);

filename = 'Test_2body_system.xlsx';
sheet = 1;
xlRange = 'E18:F27';

[v,T,vT] = xlsread(filename, sheet, xlRange);
time_step_RK4_earth_sun=v(:,1); 
distance_RK4_earth_sun=v(:,2);

filename = 'Test_2body_system.xlsx';
sheet = 1;
xlRange = 'I18:J27';

[v,T,vT] = xlsread(filename, sheet, xlRange);
time_step_VV_earth_sun=v(:,1); 
distance_VV_earth_sun=v(:,2);

%Plot data
%distance as fct of time step. Sun not moving
figure
plot(time_step_RK4_earth,distance_RK4_earth,'-.b',time_step_VV_earth,distance_VV_earth,'m','LineWidth',2)
legend('Runge-Kutta 4','Velocity-Verlet','Location','southwest')

xlabel('Length of Time Steps (days)','FontSize',12)
ylabel('Final Distance (AU)','FontSize',12)

figure
plot(time_step_RK4_earth_sun,distance_RK4_earth_sun,'-.b',time_step_VV_earth_sun,distance_VV_earth_sun,'m','LineWidth',2)
legend('Runge-Kutta 4','Velocity-Verlet','Location','southwest')

xlabel('Length of Time Steps (days)','FontSize',12)
ylabel('Final Distance (AU)','FontSize',12)


