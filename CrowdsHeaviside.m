function [time, xsols, vsols, Mass, l] = CrowdsHeaviside(t);
tic
n=50;
r = 1/4;
mcase=1;
Case=1;
% masses
if mcase == 1
   Mass = ones(n,1);
elseif mcase == 2
   if mod(n,2)==0
       Mass = ones(n,1);
       Mass(n/2)=2;
   else
       Mass = ones(n,1);
       Mass((n+1)/2)=2;
   end
elseif mcase == 3
   if mod(n,2)==0
       Mass = ones(n,1);
       Mass(n/2)=1/2;
   else
       Mass = ones(n,1);
       Mass((n+1)/2)=1/2;
   end
    
elseif mcase == 4
    if n==2;
        Mass(1)=2;
        Mass(n)=1;
    else
    Mass(2:n-1) = ones(n-2,1);
    Mass(1)=2;
    Mass(n)=2;
    end
elseif mcase == 5
    if n==2;
        Mass(1)=1/2;
        Mass(n)=1;
    else
    Mass(2:n-1) = ones(n-2,1);
    Mass(1)=1/2;
    Mass(n)=1/2;
    end
elseif mcase == 6
    % random distribution
    Mass = 1/2+3/2*(rand(n,1));
end
%Mass = ones(n,1);
%Mass = ((1/2)+(3/2)*rand(n,1));

% set length
in = r*0;
l = 2*r*n-in;
% set time interval
time_int=[0,t];
% set initial conditions
if Case == 1
    % Bell Curve
    dVec = zeros(1,n/2);
    mid = (r/2);
    for i = 1:(n/2)
        dVec(i) = (r/2)*(1-i/n);
    end
    X0 = zeros(1,n);
    X0(n/2) = l/2-mid/2;
    for i = (n/2-1):-1:1
        X0(i) = X0(i+1)+dVec(i+1)-2*r;
    end
    X0(n/2+1) = l/2+mid/2;
    for i = (n/2+2):1:n
        X0(i) = 2*r + X0(i-1) - dVec(n+1-i);
    end
    X0 = X0';
elseif Case == 2
    % Wavy
    del0 = zeros(1,n+1);
    del0(1) = r/2;
    del0(n+1) = r/2;
    del0(3:2:n-1) = 2*r-.1;
    del0(2:2:n+1) = 0;
    
    X0 = zeros(1,n);
    X0(1) = r-del0(1);
    for i = 2:n
        X0(i) = X0(i-1)+2*r-del0(i);
    end
    X0 = X0';
    l = X0(n)+r-del0(n+1)-r;
elseif Case == 3
    % top person moved down
    X0 = [(r:2*r:2*r*(n-1)),2*r*n-3/2*r]';
    
elseif Case == 4
    % top person moved down and bottom
    % person moved up both by 2*r
    X0 = [3/2*r,(3*r:2*r:2*r*(n-1)),2*r*n-3/2*r]';
    
elseif Case == 5
    % Uniform distribution. No movement
    % expected
    X0 = (r:2*r:2*r*n)';
    
elseif Case == 6
    % random distribution
    X0 = l*sort(rand(n,1));
elseif Case == 7
    % slack
    r=1/4;
l=2*(2*r);
Del(2)=r/8;
X0(1)=l/2+Del(2)/2-r;
X0(2)=l/2-Del(2)/2+r;
X0 = X0';
elseif Case == 8
%no slack

r=1/4;
Del(2)=r/8;
l=2*2*r-Del(2);
X0(1)=l/2+Del(2)/2-r;
X0(2)=l/2-Del(2)/2+r;
X0 = X0';
elseif Case == 9
%slack

r=1/4;
l=2*(2*r);
Del(2)=r/8;
X0(1)=r/2;
X0(2)=r/2+2*r-Del(2);
X0 = X0';
elseif Case == 10

%no slack

r=1/4;
Del(2)=r/8;
l=2*(2*r)-Del(2);
X0(1)=r/2;
X0(2)=r/2+2*r-Del(2);
X0 = X0';
end
m = 1;
r0 = 2*r;
M = m/2;
%REVISAR LO DEL NUMERO DE ENTRADAS
if nargin < 4
    M = m;
    
elseif M == 0
    M = 10^-12;
end

d0 = r0*(1-m/M);

mp = 1;
Mp = mp*2;
%revisar lo del numero de entradas
if nargin < 4
    M = m;
    
elseif M == 0
    M = 10^12;
end
%este dp es correcto?
dp = r0*(1-mp/Mp);

V0 = zeros(n,1);

init_cond = [X0;V0];

% ode solver
options = odeset('RelTol',1e-3);
[time, sols]=ode45(@eqn,time_int,init_cond,options);
xsols = sols(:,(1:n));
vsols = sols(:,(n+1:2*n));
% plot solutions
figure;
hold on
plot(time,xsols);
hold off
axis([time_int,0,l]);
xlabel('Time [s]','FontSize',15,'FontName','Arial')
ylabel('Location [m]','FontSize',15,'FontName','Arial')
str = sprintf('%i people of radius %i',n,r);
title({str},'FontSize',20,'FontName','Arial');

% ode function
    function dwdt = eqn(~,w)
        
        % set w-vector to differentiate
        X = zeros(1,n);
        V = zeros(1,n);
        X(1:n) = w(1:n);
        V(1:n) = w(n+1:2*n);
        
        % make vector of deltas 
        Delh = zeros(1,n+1);
        Delh(1) = (r-X(1));
        Delh(2:n) = X(1:n-1)-X(2:n)+2*r;
        Delh(n+1) = (X(n)-l+r);
        Del = Delh.*(Delh > 0);
        % compute forces
        FN = zeros(1,n+1);
        FP = zeros(1,n+1);
        k=80*1/10;
        P=5.645;
        %Delta0=0.001;
        Delta0=-0.5;
        %NW1
        FN(1) = k*tan(2*pi*Del(1));
        %PW1
        FP(1) = 0;
        % middle people
        %normal
        FN(2:n) = k*tan(pi*Del(2:n));
        %push
        
        %FP(2:n)= (P/2)*((2/pi)*atan(100*Delh(2:n))+1).*((sign(V(1:n-1)-V(2:n))).^2).*((sign(V(1:n-1)-V(2:n))+1)/2)...
        %    +(P/2)*(2/pi)*atan(100*(Delh(2:n)-Delta0)+1).*((sign(V(1:n-1)-V(2:n))).^2).*((sign(V(1:n-1)-V(2:n))-1)/2).^2 ...
        %    +(P/2)*(2/pi)*atan(100*(Delh(2:n)-Delta0)+1).*((sign(V(1:n-1)-V(2:n))-1)/2).*((sign(V(1:n-1)-V(2:n))+1)/2);
        %                FP(2:n)= (P/2)*((2/pi)*atan(100*Delh(2:n))+1).*((V(1:n-1)-V(2:n))>0)...
        %    +(P/2)*(2/pi)*atan(100*(Delh(2:n)-Delta0)+1).*(sign(V(1:n-1)-V(2:n))<=0);
        FP(2:n)= (P/20)*((V(1:n-1)-V(2:n))>0)...
            +(P/20)*(sign(V(1:n-1)-V(2:n))<=0).*((Delh(2:n)-Delta0)>0);
        %FP(2:n)= (P/20)*((Del(2:n))>=0);
        
        %NnW
        FN(n+1) = k*tan(2*pi*Del(n+1));
        %PnW
        FP(n+1) = 0;
        
        
        
        Masss(1:n)=Mass(1:n);
        %construct ode
        dwdt = zeros(2*n,1);
        dwdt(1:n) = V(1:n);
        dwdt(n+1:2*n) = ((FP(1:n)-FP(2:n+1))+FN(1:n)-FN(2:n+1))./Masss(1:n);
    end
toc
end