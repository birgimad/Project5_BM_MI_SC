close all
clear all
clc

%Import data and define closed-form solution
%Data for ordered spins (up = 1) and T = 1.0
filename = 'Project5_Different_N.xlsx';
sheet = 1;
xlRange = 'J3:J89';
[v,T,vT] = xlsread(filename, sheet, xlRange);
distance=v(:,1);
N = 100;
n0 = 1*N^2;
r0 = 2*N^(-1/3);

%x = 5:0.01:30;
%n = n0./(1+(x./r0).^4);

%Plot histograms
figure
xbins_1 = 0:0.5:30;
[f_1,x_1] = hist(distance,xbins_1);
dx_1 = diff(x_1(1:2))
bar(x_1,f_1/(4*pi*f_1.^2*0.5),'b')    %/(4*pi*x_1.^2*0.5)
%hold on 
%plot(x,n,'r','LineWidth',2)
legend('N = 30','Location','northwest')
xlabel('Distance from cluster center (lightyears)', 'fontsize',14) 
ylabel('probability','fontsize',14)




