% tipping diagram for the plant herbivore model (O'Keeffe and Wieczorek 2020) 

warning off
clc
clear
set(0,'defaulttextInterpreter','latex');

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
EPS = 1e-1;

init = [10.7008; 7.4913];

rho_0 = 0.5;
rho_fold = 0.6966;
drho_fold = rho_fold - rho_0;

res = 100;

% dnu_scan = linspace(0,2,res);
% drho_scan = linspace(0,0.5,res);
r_th = logspace(-3,0.2,res);
tipping_mat = zeros(res,res);

%% thresholds
dnu_ram = 0.066;
dnu_ret = 0.066;

for ind_r = 1: length(r_th)
    r  = r_th(ind_r);
    tspan(1) = min(-10/r,-1000);
    tspan(2) = max (50/r, 5000);
    fun_ram = @(input) does_it_tip_ram(r,input,drho_fold,rho_0,tspan,init, EPS,opts);
    fun_ret = @(input) does_it_tip_ret(r,input,drho_fold,rho_0,tspan,init, EPS,opts);
    dnu_ram = fzero(fun_ram,dnu_ram);
    dnu_ret = fzero(fun_ret,dnu_ret);
    dnu_ram_th(ind_r) = dnu_ram;
    dnu_ret_th(ind_r) = dnu_ret;
    disp(ind_r)
end
%% BI point 
dnu_ram = 0.066;
r = inf;
fun_BI = @(input) does_it_tip_ram(r,input,drho_fold,rho_0,tspan,init, EPS,opts);
dnu_BI = fzero(fun_BI,dnu_ram);
%% plotting
figure(1);
% cla
hold on
% plott = pcolor(drho_scan,r_scan,tipping_mat);
% plott.LineStyle = 'none';
% plott.FaceAlpha = 0.3;

plot(dnu_ram_th,r_th,'--b','LineWidth',1.5)
plot(dnu_ret_th,r_th,'-b','LineWidth',1.5)

plot([dnu_BI dnu_BI], [r_th(1) r_th(end)],'--k') 
plot([1 1], [r_th(1) r_th(end)],'--k') 

colormap([1 0 0; 0 1 0; 1 1 1]);
set(gca, 'YScale', 'log')
set(gcf, 'renderer', 'painters')
axis([0 1.4 1e-3 1])
%% testing 
% opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
% drho = (0.125+0.120241)/2;
% r    = 0.01;
% rho  = @(tt)((rho_0 + drho * sech(r*tt))*(tt<0) + (rho_0 + drho)*(tt>=0));
% odefun = @(tt,x)ecoODE(tt,x,rho);
% [tt,var]   = ode45(odefun,[-1000,1000],init,opts);
% figure(3);
% hold on
% plot(tt,var(:,1),...
%     tt,var(:,2))
% % figure(5); hold on
% % fplot(rho,[-100 1000])
%% fuction 
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

% tipping threshold functions 
% ramp shifr 
function output = does_it_tip_ram(r,dnu,drho,rho_0,tspan,init, EPS,opts)
% output -1 no tipping, or 1 tipping

nu_ram   = @(tt)(dnu*sech(r*tt))*(tt<0)+(dnu)*(tt>=0);
rho_ram  = @(tt)rho_0 + drho*nu_ram(tt);

odefun_ram = @(tt,x)ecoODE(tt,x,rho_ram);
[~,var_ram]   = ode45(odefun_ram,tspan,init,opts);
if (var_ram(end,2) < EPS)
    output = -1;
else
    output = 1;
end
end

% return shift
function output = does_it_tip_ret(r,dnu,drho,rho_0,tspan,init, EPS,opts)
% output -1 no tipping, or 1 tipping

nu_ret   = @(tt)(dnu*sech(r*tt));
rho_ret  = @(tt)rho_0 + drho*nu_ret(tt);

odefun_ret = @(tt,x)ecoODE(tt,x,rho_ret);
[~,var_ret]   = ode45(odefun_ret,tspan,init,opts);
if (var_ret(end,2) < EPS)
    output = -1;
else
    output = 1;
end
end