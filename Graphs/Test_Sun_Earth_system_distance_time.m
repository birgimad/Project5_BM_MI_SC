close all
clear all
clc

filename = 'Test_2body_system.xlsx';
sheet = 2;
xlRange = 'D3:F12';

[v,T,vT] = xlsread(filename, sheet, xlRange);
RK4_time_years=v(:,1); 
RK4_time_days=v(:,2);
RK4_distance_AU=v(:,3);

filename = 'Test_2body_system.xlsx';
sheet = 2;
xlRange = 'H3:J21';

[v,T,vT] = xlsread(filename, sheet, xlRange);
VV_time_years=v(:,1); 
VV_time_days=v(:,2);
VV_distance_AU=v(:,3);

%Plot data
%distance as fct of time step. Sun not moving
figure
plot(RK4_time_years,RK4_distance_AU,'-.b','LineWidth',2)
legend('Runge-Kutta 4','Location','southwest')

xlabel('Time (years)','FontSize',12)
ylabel('Final Distance (AU)','FontSize',12)

figure
plot(VV_time_years,VV_distance_AU,'m','LineWidth',2)
legend('Velocity-Verlet','Location','southwest')

xlabel('Time (years)','FontSize',12)
ylabel('Final Distance (AU)','FontSize',12)




