clear all
close all

%for i=1:6
i=1;
n=50;
t=100;
M=1;
Mp=1;
%case 1 to 5 if n=50, case 7 to 10 if n=2
Case=1;
mcase=i;
[time, xsols, vsols, Mass, l] = Crowds2(n,t,M,Mp,Case,mcase);
filename = [ 'simulacionpruebaHeaviside' num2str(i) '.mat' ];
save(filename);
clear all
close all
%end