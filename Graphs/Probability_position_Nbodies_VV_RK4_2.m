close all
clear all
clc

%Import data and define closed-form solution
%Data for ordered spins (up = 1) and T = 1.0
filename = 'Energy_and_equilibrium_study.xlsx';
sheet = 5;
xlRange = 'A4:B103';
[v,T,vT] = xlsread(filename, sheet, xlRange);
pos_initial_1=v(:,1);
pos_final_1=v(:,2);

%Plot histograms
figure
xbins_1_final = 0:0.5:1000;
[f_1_final,x_1_final] = hist(pos_final_1,xbins_1_final);
dx_1_final = diff(x_1_final(1:2));
bar(x_1_final,f_1_final/(2*sum(f_1_final*dx_1_final)),'b')
legend('Final distance, t = 4\tau_{crunch}')
xlabel('Distance from cluster center (lightyears)', 'fontsize',14) 
ylabel('probability','fontsize',14)


%Checking that procent = 1

k1_final=2*sum(f_1_final*dx_1_final);
procent1_final = sum(f_1_final/k1_final)




