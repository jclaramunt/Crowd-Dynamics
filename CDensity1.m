function [P] = CDensity1(time, xsols,l);
tic
%x step
s = 1/10;
%t step
S = 10;
% plot density
r = 1/4;
de = 1/4;
ep = 1/4;
lt = length(time);
n = length(xsols(1,:));
X = (0:s:2*r*n);
lx = length(X);
%pc = sqrt(6)*ones(1,lx);
%figure 1;
i=1:S:lt-1;
nlt = length(i);
P = zeros(nlt,lx);
k = 1;
for i=1:S:lt-1
    for j = 1:lx
        c2 = (xsols(i,:)+ep-X(j)+de).*(((xsols(i,:)-ep)<(X(j)-de)).*((xsols(i,:)+ep)>(X(j)-de)));
        c3 = (2*ep)*sum(((xsols(i,:)-ep)>(X(j)-de)).*((xsols(i,:)+ep)<(X(j)-de)));
        c4 = (-xsols(i,:)+ep+X(j)+de).*(((xsols(i,:)-ep)<(X(j)+de)).*((xsols(i,:)+ep)>(X(j)+de)));
        P(k,j) = (1/(2*de))^1*(sum(c2)+sum(c3)+sum(c4));      
    end 
    k = k + 1;
%     plot(X,P(i,:))
%     %hold on
%     %plot(X,pc)
%     %hold off
%     axis([-.1 2*r*n+.1 -.1 4])
%     str = sprintf('Density at t = %i [s]',time(i));
%     title({str},'FontSize',20,'FontName','Arial');
%     xlabel('Location [meters]','FontSize',15,'FontName','Arial')
%     ylabel('Density [people per meter]','FontSize',15,'FontName','Arial')
%     dt = time(i+1)-time(i);
%     pause(10^-250/dt)   
%     F(i) = getframe(gcf);
   
end
toc
  contour(time(1:S:lt-1),X,P')
    axis([0,max(time),0,l+0.75]);
xlabel('Time [s]','FontSize',12,'FontName','Arial')
ylabel('Location [m]','FontSize',12,'FontName','Arial')
str = sprintf('Density Plot of Dynamical Crowd System');
title({str},'FontSize',15,'FontName','Arial');
MP=max(P,[],2);
avep=mean(MP);
MAXP=max(MP);
% figure;
% hold on
% plot([1:length(MP)],MP);
% hold off
% axis([1,length(MP),0,MAXP+0.5]);
% str = sprintf('Maximum density: %i',MAXP');
% str2 = sprintf('Average Maximum density: %i',avep');
% %COMPUTATIONAL STEPS DOESNT SOUND WELL
% xlabel('Computational steps','FontSize',15,'FontName','Arial');
% ylabel('Maximum density','FontSize',15,'FontName','Arial');
% annotation('textbox', [0.05,0.00,0.1,0.1],...
%            'String',{str2},'FontSize',11);
% annotation('textbox', [0.05,0.88,0.1,0.1],...
%            'String',{str},'FontSize',11);
end