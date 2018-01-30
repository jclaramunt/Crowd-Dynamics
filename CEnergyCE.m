function [E,TKE] = CEnergyCE(xsols,vsols,time,Mass,l);
tic
%--------------------------------------------------------%
% initial parameters
r = 1/4;
in = r*0;
[lt, n] = size(xsols);
t = time(lt);
time_int=[0,t];
%l = 2*r*n-in;
k=80*1/10;
M =(Mass').^-1;
%--------------------------------------------------------%
% compute Kinetic Energies
%--------------------------------------------------------%
P=5.645;
KEi =(vsols.^2)/2;
% divide by masses
KE = bsxfun(@rdivide,KEi',M);
TKE = sum(KE)';
%--------------------------------------------------------%
% compute Potential Energies
%--------------------------------------------------------%
% get deltas
Del = zeros(lt,n+1);
for i = 1:lt
    Del(i,1) = (r-xsols(i,1));
    for j = 2:n;
        Del(i,j) = (xsols(i,j-1)+r)-(xsols(i,j)-r);
    end
    Del(i,n+1) = (xsols(i,n)-l+r);    
end

%con la nueva definicion de los delta esto cambia
Del = Del.*(Del > 0);



% plug into potential function
%PE = -(4*r*k/(2*pi))*log(abs(cos(pi*Del/(4*r))));
%PE(:,1) = -(2*r*k/(2*pi))*log(abs(cos((2*pi)*Del(:,1)/(2*r))));
%PE(:,n+1) = -(2*r*k/(2*pi))*log(abs(cos((2*pi)*Del(:,n+1)/(2*r))));
%PE = -(k/pi)*log(abs(cos(pi*Del)));
%PE(:,1) = -(k/(2*pi))*log(abs(cos((2*pi)*Del(:,1))));
%PE(:,n+1) = -(k/(2*pi))*log(abs(cos((2*pi)*Del(:,n+1))));
PE = -(4*k*r/pi)*log(abs(cos(pi*Del/(4*r))));
PE(:,1) = -(2*k*r/(pi))*log(abs(cos((pi)*Del(:,1)/(2*r))));
PE(:,n+1) = -(2*k*r/(pi))*log(abs(cos((pi)*Del(:,n+1)/(2*r))));
% sum
TPE_init = sum(PE')';
%Pmx = max(TPE_init);
Pmn = min(TPE_init);
% fix offset
TPE = TPE_init - Pmn;
%--------------------------------------------------------%
% Total Energy
E = TKE+TPE;
mx = max(E);
mn = min(E);
ave = (mx+mn)/2;
%--------------------------------------------------------%
% plot solutions
figure
subplot(2,1,1);
plot(time,E,'m','LineWidth',2.0);

% M1 = 'Total';
% legend({M1},'FontSize',16,...
%     'Location','southoutside','Orientation','horizontal');
str = sprintf('Average Energy %i (J*)',ave);
annotation('textbox', [0.0,0.04,0.1,0.1],...
           'String',{str},'FontSize',11);
axis([time_int,0-.01,1.1*mx]);
xlabel('Time (s)','FontSize',15,'FontName','Arial')
ylabel('Energy (J*)','FontSize',15,'FontName','Arial')
title('Energy Growth of Dynamic Crowd System','FontSize',20,'FontName','Arial');


subplot(2,1,2)
hold on
plot(time,E,'m','LineWidth',2.0);
plot(time,TKE,'b','LineWidth',0.2);
plot(time,TPE,'r','LineWidth',0.2);
M1 = 'Total';
M2 = 'KE';
M3 = 'PE';
legend({M1, M2, M3},'FontSize',16,...
    'Location','southoutside','Orientation','horizontal');
% str = sprintf('Average Energy %i [J]',ave);
% annotation('textbox', [0.0,0.04,0.1,0.1],...
%            'String',{str},'FontSize',11);
axis([time_int,0-.01,1.1*mx]);
xlabel('Time (s)','FontSize',15,'FontName','Arial')
ylabel('Energy (J*)','FontSize',15,'FontName','Arial')
title('Energy trend of Dynamic Crowd System','FontSize',20,'FontName','Arial');
hold off
%--------------------------------------------------------%
toc
end

