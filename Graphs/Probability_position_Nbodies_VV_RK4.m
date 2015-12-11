close all
clear all
clc

%Import data and define closed-form solution
%Data for ordered spins (up = 1) and T = 1.0
filename = 'Time_step_length_N_bodies_RK4.xlsx';
sheet = 9;
xlRange = 'A3:B102';
[v,T,vT] = xlsread(filename, sheet, xlRange);
pos_initial_1=v(:,1);
pos_final_1=v(:,2);
filename = 'Time_step_length_N_bodies_VV.xlsx';
sheet = 5;
xlRange = 'A2:B101';
[v,T,vT] = xlsread(filename, sheet, xlRange);
pos_initial_2=v(:,1);
pos_final_2=v(:,2);
filename = 'Time_step_length_N_bodies_VV.xlsx';
sheet = 3;
xlRange = 'A2:B101';
[v,T,vT] = xlsread(filename, sheet, xlRange);
pos_initial_3=v(:,1);
pos_final_3=v(:,2);

%Plot histograms
figure
xbins_1 = 0:0.5:20;
[f_1,x_1] = hist(pos_initial_1,xbins_1);
dx_1 = diff(x_1(1:2));
bar(x_1,f_1/(2*sum(f_1*dx_1)),'r')
hold on 
xbins_1_final = 0:0.5:20;
[f_1_final,x_1_final] = hist(pos_final_1,xbins_1_final);
dx_1_final = diff(x_1_final(1:2));
bar(x_1_final,f_1_final/(2*sum(f_1_final*dx_1_final)),'b')
legend('Initial distance','Final distance, dt = 10^3 yr','Location','northwest')
xlabel('Distance from cluster center (lightyears)', 'fontsize',14) 
ylabel('probability','fontsize',14)

figure
xbins_2 = 0:0.5:20;
[f_2,x_2] = hist(pos_initial_2,xbins_2);
dx_2 = diff(x_2(1:2));
bar(x_2,f_2/(2*sum(f_2*dx_2)),'r')
hold on 
xbins_2_final = 0:0.5:20;
[f_2_final,x_2_final] = hist(pos_final_2,xbins_2_final);
dx_2_final = diff(x_2_final(1:2));
bar(x_2_final,f_2_final/(2*sum(f_2_final*dx_2_final)),'b')
legend('Initial distance','Final distance, dt = 100 yr','Location','northwest')
xlabel('Distance from cluster center (lightyears)', 'fontsize',14) 
ylabel('probability','fontsize',14)

figure
xbins_3 = 0:0.5:20;
[f_3,x_3] = hist(pos_initial_3,xbins_3);
dx_3 = diff(x_3(1:2));
bar(x_3,f_3/(2*sum(f_3*dx_3)),'r')
hold on 
xbins_3_final = 0:0.5:20;
[f_3_final,x_3_final] = hist(pos_final_3,xbins_3_final);
dx_3_final = diff(x_3_final(1:2));
bar(x_3_final,f_3_final/(2*sum(f_3_final*dx_3_final)),'b')
legend('Initial distance','Final distance, dt = 10^4 yr','Location','northwest')
xlabel('Distance from cluster center (lightyears)', 'fontsize',14) 
ylabel('probability','fontsize',14)

%Checking that procent = 1

k1=2*sum(f_1*dx_1);
procent1 = sum(f_1/k1)
k2=2*sum(f_2*dx_2);
procent2 = sum(f_2/k2)
k3=2*sum(f_3*dx_3);
procent3 = sum(f_3/k3)

k1_final=2*sum(f_1_final*dx_1_final);
procent1_final = sum(f_1_final/k1_final)
k2_final=2*sum(f_2_final*dx_2_final);
procent2_final = sum(f_2_final/k2_final)
k3_final=2*sum(f_3_final*dx_3_final);
procent3_final = sum(f_3_final/k3_final)



