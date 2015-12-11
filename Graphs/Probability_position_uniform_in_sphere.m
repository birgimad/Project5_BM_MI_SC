close all
clear all
clc

%Import data and define closed-form solution
%Data for ordered spins (up = 1) and T = 1.0
filename = 'Project5_random_mass_check.xlsx';
sheet = 3;
xlRange = 'A1:B100000';

[v,T,vT] = xlsread(filename, sheet, xlRange);
data_r=v(:,1); 
%data_final=v(:,2); 
%data_phi=v(:,5);
%data_theta=v(:,7);

%Plot histograms
figure

xbins_r = -19.5:0.5:19.5;
onevector = ones(size(xbins_r));

[f_r,x_r] = hist(data_r,xbins_r);

dx_r = diff(x_r(1:2));
bar(x_r,(f_r./(pi*(400-(xbins_r).^2))),'r')

xlabel('Distance in x-direction from cluster center (ly)', 'fontsize',14) % x-axis label
ylabel('Density in cross-sectional area (ly ^{-2})','fontsize',14) % y-axis label




