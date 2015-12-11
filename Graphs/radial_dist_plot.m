close all
clear all
clc

filename = 'Different_N.xlsx';
sheet = 1;
xlRange = 'A4:B7';

[v,T,vT] = xlsread(filename, sheet, xlRange);
N=v(:,1); 
R=v(:,2);

n0 = .^2;
r0 = .^(-1/3);

x = 12:0.01:15;
%n = n0/(1+(x/r0).^4);

%Plot data
figure
plot(R,N,'bx','LineWidth',2)
hold on 
plot(x,n,'r','LineWidth',2)
xlabel('Mean distance to cluster center (ly)','FontSize',12)
ylabel('Number of particles','FontSize',12)





