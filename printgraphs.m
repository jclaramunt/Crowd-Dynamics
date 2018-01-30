r=0.25;
t=100;
%t=100;
time_int=[0,t];
figure;
hold on
plot(time,xsols);
hold off
axis([time_int,0,l+0.5]);
xlabel('Time (s)','FontSize',15,'FontName','Arial')
ylabel('Location (m)','FontSize',15,'FontName','Arial')
str = sprintf('%i people of radius %i',n,r);
title({str},'FontSize',20,'FontName','Arial');