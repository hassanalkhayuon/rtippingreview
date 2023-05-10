% basin instability analysis for the AMOC 3 box model system. 

warning off
clc
clear
addpath(...
    'C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

opts = odeset('RelTol',1e-5,'AbsTol',1e-10);

tspan = [0,3000];

H_start     = 0.25;
G_start     = 0.36;

parStart    = [H_start,G_start];


initStabEq  = [0.035435;0.034912];
gg          = @(var)AMOC_autoODE(var,parStart);
[stabEqStart,~,conv] = NR(gg,initStabEq);

[~,var] = ode45(@(t,var)gg(var),tspan,stabEqStart,opts);

Hres  = 500;
Gres  = 500;

H_scan = linspace(0,0.42,Hres);
G_scan  = linspace(-4,2,Gres);

scanMatrix = NaN(Hres,Gres);

% Define Fe1, Fe2_upp, Fe2_low and H curves. 
load('data/2parBD_AMOC.mat')

n = 1;
for ind_hopf = 1:length(Hopf)
    if mod(ind_hopf,10) == 0
        Hopfn(n,:) = Hopf(ind_hopf,:);
        n = n+1;
    end
end

gamma_F     = @(HH)interp1(Fold(:,1),Fold(:,2),HH);
gamma_H     = @(QQ)interp1(Hopfn(:,1),Hopfn(:,2),QQ);


%% scan the 2 parameter plane
for ind_h = 1:Hres
    for ind_g = 1:Gres
        
        H     = H_scan(ind_h);
        gamma = G_scan(ind_g);
        par = [H,gamma];
        
        odefun = @(t,var)AMOC_autoODE(var,par);
        [~,var] = ode45(odefun,tspan,stabEqStart,opts);
        
        if ~isnan(gamma_H(H))
            parSpaceCondition = gamma >= max(gamma_H(H),gamma_F(H));
        else
            parSpaceCondition = gamma >= gamma_F(H);
        end
        tippingCond = var(end,1)<-0.1;
        
        if(parSpaceCondition  && tippingCond)
            scanMatrix(ind_g,ind_h) = 1;
        end
         disp([ind_g,ind_h]);
    end
end
%% calculate the threshold
 
% % H = 0
% H_the = linspace(0.2735,0.4032,100);
% ginit = -4;

% % H = 0.2
% H_the = linspace(0.2752,0.4141,100);
% ginit = -4;

% H = 0.25
H_the = linspace(0.2752,0.4183,100);
ginit = -4;

% % H = 0.3;
% H_the = linspace(0.2769,0.415,100);
% ginit = -3.988;


g_the = NaN(size(H_the));

for ind_H = 1:length(H_the)
    H = H_the(ind_H);
    therFun = @(gg)ThreFun(H,gg,stabEqStart);
    ginit = fzero(therFun,ginit);
    g_the(ind_H) = ginit;
    disp([ind_H,100])
end 

%% plotting 
grayscal = 1:-0.1:0;
map = [...
    1.0,1.0,1.0;...
    0.8,0.8,0.8;...
    0.6,0.6,0.6];
figure; hold on
PCOLOR = pcolor(H_scan,G_scan,scanMatrix);
PCOLOR.FaceAlpha = 1;
PCOLOR.LineStyle = 'none';
colormap(map)
caxis([0 2]);
hold on
fplot(gamma_H)
fplot(gamma_F)
%% testing 
 
H     = 0.3619;
gamma = -0.621;
par = [H,gamma];

odefun = @(t,var)AMOC_autoODE(var,par);
[t,var] = ode45(odefun,[0 2000],stabEqStart,opts);
figure; hold on
plot(...
     t,var(:,1),'k-',...
     t,var(:,2),'b-'...
     )
  
%% functions 
% H = 0.2769;
% ginit = -3.988;
% ThreFun(H,ginit,stabEqStart)
 function [output] = ThreFun(H,gamma,stabEqStart)

opts = odeset('RelTol',1e-5,'AbsTol',1e-10);
tspan = [0,3000];

par = [H,gamma];


odefun = @(t,var)AMOC_autoODE(var,par);
[~,var] = ode45(odefun,tspan,stabEqStart,opts);

% 
% gg          = @(var)AMOC_autoODE(var,par);
% [refPoint,~,~] = NR(gg,stabEqStart);

if var(end,1)<-0.1
    output = -1;
else
    output = 1;
end
end

% function [output] = Bout_ThreFun(q,p,stabEqStart)
% 
% % find the initial condition
% opts = odeset('RelTol',1e-5,'AbsTol',1e-10,'Events', @myEvent);
% tspan = [0,100];
% par = [q,p];
% odefun = @(t,var)PowSys(var,par);
% [~,var] = ode45(odefun,tspan,stabEqStart,opts);
% 
% if (norm(var(end,4))> 1e-3)
%     output = -1;
% else
%     output = 1;
% end
% end

% the ode fuction of AMOC model
function [ z ] = AMOC_autoODE(var,par)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

S0 = 0.035; % dimensionless

alpha = 0.12; %kg m^-3 C^-1
beta = 790.0;%kg m^-3 
Y = 100*3.15e7; %sec/year

H = par(1);
gamma = par(2);


% 2 x CO2
VN = 0.3683e17;%m^3
VT = 0.5418e17; %m^3
VS = 0.6097e17; %m^3
VIP = 1.4860e17; %m^3
VB = 9.9250e17; %m^3

Ts = 7.919; %C
T0 = 3.87; %C

C = 4.4735e16; %m^3 (total sum of V*S)

lambda = 1.62e7; %m^6kg^-1s^-1
mu = 22e-8; %deg s m^-3

KN = 1.762e6; %m^3s^-1
KS = 1.872e6; %m^3s^-1

FN = (0.486+0.1311.*H).*1e6;
FT = (-0.997+0.6961.*H).*1e6; %m^3s^-1
%====================================================

SN = var(1,:);
ST = var(2,:);
SS = (0.034427-S0).*100;
%SIP = (0.034668-S0).*100;
SB = (0.034538-S0).*100;
SIP = 100.*(C-(VN.*SN+VT.*ST+VS.*SS+VB.*SB)./100-S0.*(VB+VN+VT+VIP+VS))./VIP;

q = lambda.*(alpha.*(Ts-T0)+beta.*(var(1,:)./100-SS./100))/(1+lambda*alpha*mu);
aq = abs(q);


z1p = (Y./VN).*(q.*(ST./100-SN./100)+KN.*(ST./100-SN./100)-FN.*S0);
   
z2p = (Y./VT).*(q.*(gamma.*SS./100+(1-gamma).*SIP./100-ST./100)+...
    KS.*(SS./100-ST./100)+KN.*(SN./100-ST./100)-FT.*S0);
    

z1n = (Y./VN).*(aq.*(SB./100-SN./100)+KN.*(ST./100-SN./100)-FN.*S0);
    
z2n = (Y./VT).*(aq.*(SN./100-ST./100)+KS.*(SS./100-ST./100)...
    +KN.*(SN./100-ST./100)-FT.*S0);


z(1,:) = z1p.*(q>=0)+z1n.*(q<0);
z(2,:) = z2p.*(q>=0)+z2n.*(q<0);


end

function [value, isterminal, direction] = myEvent(t,y)
TOL = 0;
value      = (y(4)<=TOL);
isterminal = 1;   % Stop the integration
direction  = 0;   % approch zero from either diractions? 
end

function [value, isterminal, direction] = threFun(q,p)
TOL = 0;
value      = (y(4)<=TOL);
isterminal = 1;   % Stop the integration
direction  = 0;   % approch zero from either diractions? 
end
