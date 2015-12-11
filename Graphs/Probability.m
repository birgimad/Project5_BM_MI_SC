close all
clear all
clc

%Import data and define closed-form solution
%Data for ordered spins (up = 1) and T = 1.0
filename = 'Project5_random_mass_check.xlsx';
sheet = 1;
xlRange = 'A1:A100000';

[v,T,vT] = xlsread(filename, sheet, xlRange);
data_random_mass=v(:,1); 

y = [7:.1:13];
norm = normpdf(y,10,1);

%Plot histograms
figure

xbins = 7:0.5:13;

[f,x] = hist(data_random_mass/10000,xbins)
dx = diff(x(1:2));
bar(x,f/sum(f*dx))

hold on

plot(y,norm,'r','LineWidth',2)

legend('Generated masses','Gaussian dist.')

xlabel('Mass in unit of solar masses', 'fontsize',14) % x-axis label
ylabel('Probability','fontsize',14) % y-axis label




