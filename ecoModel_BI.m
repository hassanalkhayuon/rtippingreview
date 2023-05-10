% basin instability analysis for the Plant Hibiour model 
% (O'keefe and wieczorek 2019) 

warning off
clc
clear
addpath(...
    'C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% model parameters
C =0.02;
E =0.4;
b = 0.025;
Cmax =1;
a =10;
bc = 0.025;

%m_tr = @(rr)( (E*Cmax*exp((-b-bc)*rr/C))/(1+(a*C/rr)^2));
% fplot(m_tr,[0,1])

opts = odeset('RelTol',1e-5,'AbsTol',1e-10);

tspan = [0,500];

r_start     = 0.5;
m_start     = 0.125;

parStart    = [r_start,m_start];


initStabEq  = [7.58;11.25];
% gg          = @(var)ecoAutODE(var,parStart);
% [stabEqStart,~,conv] = NR(gg,initStabEq);
stabEqStart = initStabEq;
% [~,var] = ode45(@(t,var)gg(var),tspan,stabEqStart,opts);

rres  = 1000;
mres  = 1000;

r_scan = linspace(0,1.5,rres);
m_scan = linspace(0.05,0.15,mres);
r_tr1   = linspace(0,.2766,100);
r_tr2   = linspace(0.2766,2,100);
m_tr1   = (E.*Cmax.*exp((-b-bc).*r_tr1./C))./(1+(a.*C./r_tr1).^2);
m_tr2   = (E.*Cmax.*exp((-b-bc).*r_tr2./C))./(1+(a.*C./r_tr2).^2);

scanMatrix = NaN(rres,mres);

% Define Fe1, Fe2_upp, Fe2_low and H curves. 
load('data/2parBD_eco.mat')

% n = 1;
% for ind_hopf = 1:length(Hopf)
%     if mod(ind_hopf,10) == 0
%         Hopfn(n,:) = Hopf(ind_hopf,:);
%         n = n+1;
%     end
% end

m_F     = @(HH)interp1(Fold(:,1),Fold(:,2),HH);
m_H     = @(QQ)interp1(Hopf(:,1),Hopf(:,2),QQ);
m_h     = @(QQ)interp1(homo(:,1),homo(:,2),QQ);

%% scan the 2 parameter plane
for ind_h = 364:rres
    for ind_g = 1:mres
        
        r     = r_scan(ind_h);
        m     = m_scan(ind_g);
        par = [r,m];
        
        odefun = @(t,var)ecoAutODE(var,par);
        [~,var] = ode45(odefun,tspan,stabEqStart,opts);
        
        if ~isnan(m_H(r))
            parSpaceCondition = m <= min(m_H(r),m_F(r));
        else
            parSpaceCondition = m <= m_F(r);
        end
        
        tippingCond = norm(var(end,2))< 1e-1;
        
        if(tippingCond && parSpaceCondition)
            scanMatrix(ind_g,ind_h) = 1;
        end
         disp([ind_g,ind_h]);
    end
end
% %% calculate the threshold
%  
% % % H = 0
% % H_the = linspace(0.2735,0.4032,100);
% % ginit = -4;
% 
% % % H = 0.2
% % H_the = linspace(0.2752,0.4141,100);
% % ginit = -4;
% 
% % H = 0.25
% H_the = linspace(0.2752,0.4183,100);
% ginit = -4;
% 
% % % H = 0.3;
% % H_the = linspace(0.2769,0.415,100);
% % ginit = -3.988;
% 
% 
% g_the = NaN(size(H_the));
% 
% for ind_H = 1:length(H_the)
%     r = H_the(ind_H);
%     therFun = @(gg)ThreFun(r,gg,stabEqStart);
%     ginit = fzero(therFun,ginit);
%     g_the(ind_H) = ginit;
%     disp([ind_H,100])
% end 

%% plotting 
grayscal = 1:-0.1:0;
map = [...
    1.0,1.0,1.0;...
    0.8,0.8,0.8;...
    0.6,0.6,0.6];
figure; hold on
PCOLOR = pcolor(r_scan,m_scan,scanMatrix);
PCOLOR.FaceAlpha = 1;
PCOLOR.LineStyle = 'none';
colormap(map)
caxis([0 2]);

hold on
plot(...
    Hopf(:,1),Hopf(:,2),'-r',...
    Fold(210:end,1),Fold(210:end,2),'-b',...
    homo(:,1),homo(:,2),'--k',...
    'LineWidth',2)

plot(...
    r_tr1,m_tr1,'-','Color',[.4 .8 0],...
    'LineWidth',2)
plot(...
    r_tr2,m_tr2,'-','Color',[.4 .8 0],...
    'LineWidth',2)
plot(...
    parStart(1),parStart(2),'k.',...
    'MarkerSize',30)
plot(...
    0.2766,0.13155,'.k',...
    'MarkerSize',30) %bitchfork point
plot(...
    0.6971,0.13155,'.k',...
    'MarkerSize',30) %BT point.
plot(...
    [0.55 1.45],[0.105 0.105],'k-',...
    'LineWidth',1)
plot(...
    [0.55 0.55],[0.105 0.145],'k-',...
    'LineWidth',1)
plot(...
    [r_start 1.45],[m_start 0.145],'k-',...
    'LineWidth',1)
% aa = 0.03870;
% bb = 0.083715;
% xx = linspace(0.55,1.3,100);
% yy = aa.*xx + bb;
% plot(...
%     xx,yy,'k-',...
%     'LineWidth',1)
box on
xlabel('$\rho$')
ylabel('$m$','Rotation',0)
set(gca,'FontSize',16)
axis([0 1.503 0.08 0.15])
%% testing 
%  
% r     = 0.3619;
% m = -0.621;
% par = [r,m];
% 
% odefun = @(t,var)ecoAutODE(var,par);
% [t,var] = ode45(odefun,[0 2000],stabEqStart,opts);
% figure; hold on
% plot(...
%      t,var(:,1),'k-',...
%      t,var(:,2),'b-'...
%      )
%   
%% functions 

% ThreFun(H,ginit,stabEqStart)
%  function [output] = ThreFun(H,gamma,stabEqStart)
% 
% opts = odeset('RelTol',1e-5,'AbsTol',1e-10);
% tspan = [0,3000];
% 
% par = [H,gamma];
% 
% 
% odefun = @(t,var)ecoAutODE(var,par);
% [~,var] = ode45(odefun,tspan,stabEqStart,opts);
% 
% % 
% % gg          = @(var)ecoAutODE(var,par);
% % [refPoint,~,~] = NR(gg,stabEqStart);
% 
% if var(end,1)<-0.1
%     output = -1;
% else
%     output = 1;
% end
% end

% the ode fuction 
function [ dvar ] = ecoAutODE(var,par)
%ecoAutODE the ode function for the eco model (O'Keeffe, Wieczorek 2019)


% bifurcation parameters
r = par(1);
m = par(2);

% other parameters
C =0.02;
E =0.4;
b = 0.025;
Cmax =1;
a =10;
bc = 0.025;

P = var(1);
H = var(2);

G = Cmax*( P^2 / (P^2 + a^2) )*exp(-bc*P);
%odes

dP = r*P - C*(P^2) - H*G;
dH = (G*E*exp(-b*P) - m)*H;

dvar = [dP;dH];
end

