% basin instability analysis for the Plant Hibiour model 
% (O'keefe and wieczorek 2019) 

warning off
clc
clear
set(0,'defaulttextInterpreter','latex');

opts = odeset('RelTol',1e-5,'AbsTol',1e-10);
load('data/1parBD_eco.mat')

Cs_bran1  = @(input)interp1(Cs(:,1),Cs(:,2),input);
Cu_bran1  = @(input)interp1(Cu(:,1),Cu(:,2),input);
H_bran1   = @(input)interp1(H(:,1) , H(:,2),input);

Cs_bran2  = @(input)interp1(Cs(:,1),Cs(:,3),input);
Cu_bran2  = @(input)interp1(Cu(:,1),Cu(:,3),input);
H_bran2   = @(input)interp1(H(:,1) , H(:,3),input);

%%
tspan_slow = [-200,150];
tspan_fast = [-200,150];
tspan_crit = [-200,150];
 
rho_0 = 0.5;
rho_fold = 0.6966;
drho_fold = rho_fold - rho_0;
dnu  = 0.5;
drho = (dnu/5 + 0.5) - rho_0;



r_slow = 0.03; %(track)
r_crit = 0.0453295; %critical), 
r_fast = 0.1; %(tipping).
%

nu_slow =  @(tt)(dnu*sech(r_slow*tt))*(tt<0)+(dnu)*(tt>=0);
nu_crit =  @(tt)(dnu*sech(r_crit*tt))*(tt<0)+(dnu)*(tt>=0);
nu_fast =  @(tt)(dnu*sech(r_fast*tt))*(tt<0)+(dnu)*(tt>=0);



rho_slow   = @(tt)rho_0 + drho_fold*nu_slow(tt);
rho_crit   = @(tt)rho_0 + drho_fold*nu_crit(tt);
rho_fast   = @(tt)rho_0 + drho_fold*nu_fast(tt);

init     = [Cs_bran1(rho_0),Cs_bran2(rho_0)]; 

odefun_slow       = @(tt,x)ecoODE(tt,x,rho_slow);
[time_slow,x_slow]        = ode45(odefun_slow,tspan_slow,init,opts);

odefun_crit       = @(tt,x)ecoODE(tt,x,rho_crit);
[time_crit,x_crit]        = ode45(odefun_crit,tspan_crit,init,opts);

odefun_fast       = @(tt,x)ecoODE(tt,x,rho_fast);
[time_fast,x_fast]        = ode45(odefun_fast,tspan_fast,init,opts);


%% Coloring the other basin of attraction 

% res_x = 500;
% res_y = res_x;
% 
% x_scan      = linspace(5,25,res_x);
% y_scan      = linspace(6,11,res_y);
% 
% Rho         = @(tt)(rho_0+Delta_r);
% Cs_EQ       = [Cs_bran1(Rho(0)),Cs_bran2(Rho(0))];
% Cu_EQ       = [Cu_bran1(Rho(0)),Cu_bran2(Rho(0))];
% H_EQ        = [H_bran1(Rho(0)),H_bran2(Rho(0))];
% 
% tspan_basin = [0 1000];
% colormat    = ones(res_y,res_x);
%  
% for ind_x = 1:res_x
%      parfor ind_y = 1:res_y
%          init_basin         = [x_scan(ind_x);y_scan(ind_y)];
%          odefun_basin       = @(tt,x)ecoODE(tt,x,Rho);
%          [~,x_basin]        = ode45(odefun_basin,tspan_basin,init_basin);
%          if norm(x_basin(end,:) - Cs_EQ) < 1
%              colormat(ind_y,ind_x) = 0;
%          end
%      end
%      disp(ind_x);
%  end
load('data/basin_max_eco.mat')

%% plotting

red = [0.86,0.02,0.05];
yellow = [221,170,51]/255;
green = [0.31,0.70,0.40];
blue = [0.10,0.40,0.69];

Cs_max      = [Cs_bran1(rho_0+drho) ; Cs_bran2(rho_0+drho)];
Cu_max      = [Cu_bran1(rho_0+drho) ; Cu_bran2(rho_0+drho)];
H_max       = [H_bran1(rho_0+drho)  ; H_bran2(rho_0+drho)];

Cs_min      = [Cs_bran1(rho_0)  ;   Cs_bran2(rho_0)];
Cu_min      = [Cu_bran1(rho_0)  ;   Cu_bran2(rho_0)];
P_min       = [H_bran1(rho_0)   ;   H_bran2(rho_0)];

rho_branch = linspace(rho_0,(rho_0+drho),100);

Cs_branch = [Cs_bran1(rho_branch); Cs_bran2(rho_branch)];

[~,w_s_max1] = ode45 (@(t,x)ecoODE(t,x,@(t)rho_0+drho), ...
    [500 0],Cu_max+[0.01; 0],opts);
[~,w_s_max2] =  ode45 (@(t,x)ecoODE(t,x,@(t)rho_0+drho), ...
    [500 0],Cu_max-[0.01; 0],opts);

[~,w_s_min1] = ode45 (@(t,x)ecoODE(t,x,@(t)rho_0), ...
    [500 0],Cu_min+[0.01; 0],opts);
[~,w_s_min2] =  ode45 (@(t,x)ecoODE(t,x,@(t)rho_0), ...
    [500 0],Cu_min-[0.01; 0],opts);

figure(1); cla
hold on
grayscal = 1:-0.01:0;
map = [grayscal',grayscal', grayscal'];
PCOLOR = pcolor(x_scan,y_scan,colormat);
PCOLOR.FaceAlpha = 0.1;
PCOLOR.LineStyle = 'none';
colormap(map)

% tipping threshold W^s(Cu) at max forcing 
plot(...
    w_s_max1(:,1),w_s_max1(:,2),'-',...
    w_s_max2(:,1),w_s_max2(:,2),'-',...
    'Color',[0.7 0.7,0.7],'LineWidth',1.5)

% tipping threshold W^s(Cu) at min forcing 
plot(...
    w_s_min1(:,1),w_s_min1(:,2),'-',...
    w_s_min2(:,1),w_s_min2(:,2),'-',...
    'Color',[0.7 0.7,0.7],'LineWidth',1.5)




% critical trajectory 1
plot(x_crit(:,1),x_crit(:,2),'-',...
    'Color',blue,'LineWidth',2)

% slow trajectory 
plot(x_slow(:,1),x_slow(:,2),'-',...
    'Color',green,'LineWidth',2)

% % fast trajectory
plot(x_fast(:,1),x_fast(:,2),'-',...
    'Color',red,'LineWidth',2)

% the stable brach of the coexectance state
plot(...
    Cs_branch(1,:),Cs_branch(2,:),':k',...
    'Color','k','LineWidth',1.5);


% the stable coexestance state at max forcing 
plot(...
    Cs_max(1),Cs_max(2),'.k',...
    'MarkerSize',30)

% the stable coexestance state at min forcing
plot(...
    Cs_min(1),Cs_min(2),'.k',...
    'MarkerSize',30)

% the unstable coexestance state at max forcing
plot(...
    Cu_max(1),Cu_max(2),'ok',...
    'MarkerSize',7,'LineWidth',2)

% the stable Plant only state at min forcing
plot(...
    H_max(1),H_max(2),'.k',...
    'MarkerSize',30)

box on
% axis([0 40 0 15])
set(gca,'FontSize',16)

%% time sieries plots 

figure(10);
cla
hold on

% critical trajectory 1
plot(time_crit,x_crit(:,2),'-',...
    'Color',blue,'LineWidth',2)

% slow trajectory 
plot(time_slow,x_slow(:,2),'-',...
    'Color',green,'LineWidth',2)

% % fast trajectory
plot(time_fast,x_fast(:,2),'-',...
    'Color',red,'LineWidth',2)

figure(11);
cla
hold on

fplot (nu_crit,[-200 150],'Color',blue,'LineWidth',2)
fplot (nu_slow,[-200 150],'Color',green,'LineWidth',2)
fplot (nu_fast,[-200 150],'Color',red,'LineWidth',2)
%% the ode fuction 
function [ dvar ] = ecoODE(t,var,par)

% bifurcation parameters
rho = par(t);

g = 0.005/0.15;
c = 0.125 - g*0.5;

% Mortality rate, m, is linearly linked to plant growth rate, rho, by
m = g*rho + c;

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

dP = rho*P - C*(P^2) - H*G;
dH = (G*E*exp(-b*P) - m)*H;

dvar = [dP;dH];
end

