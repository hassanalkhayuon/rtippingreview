% tipping diagram for the AMOC model (Alkhayuon et al 2019) 
% define on and off state for dH and compare with it!
warning off
clc
clear
set(0,'defaulttextInterpreter','latex');

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

%% find off- and on-state 
load('data/1parBD_AMOC.mat')
ds(:,3) = ds(:,3) - 0.0029;

up_bran1  = @(input)interp1([us(1:end-1,1);uu(:,1)],...
    [us(1:end-1,2);uu(:,2)],input);
lo_bran1 = @(input)interp1(ds(:,1),ds(:,2),...
    input);
un_bran1 = @(input)interp1(du(:,1),du(:,2),...
    input);

up_bran2  = @(input)interp1([us(1:end-1,1);uu(:,1)],...
    [us(1:end-1,3);uu(:,3)],input);
lo_bran2 = @(input)interp1(ds(:,1),ds(:,3),...
    input);
un_bran2 = @(input)interp1(du(:,1),du(:,3),...
    input);

%%

EPS = 1e-1;
% tspan = [-1000000, 1000000];
init = [0.0324    0.1434];

H0 = 0;

res = 100;
dH_scan = linspace(0.35,0.45,res);
% dH_scan = linspace(0.412,0.4165,res);
r_scan = logspace(-6,-1,res);
% r_scan = linspace(3e-5,5e-5,res);
tipping_mat = NaN(res,res);

%% scan 
% for ind_r = 1: length(r_scan)
%     for ind_drho = 1:length(dH_scan)
%         r        = r_scan(ind_r);
%         dH       = dH_scan(ind_drho);
%         AMOC_on  = [up_bran1(dH), up_bran2(dH)]; 
%         AMOC_off = [lo_bran1(dH), lo_bran2(dH)];
%         tip_ram  = does_it_tip_ram(r,dH,H0,tspan,init,opts);
%         tip_ret  = does_it_tip_ret(r,dH,H0,tspan,init,opts);
%         tipping_mat(ind_r, ind_drho) = tip_ram + tip_ret;
%     end
%     disp([ind_drho,ind_r])
% end
%% thresholds
fzeroopts = optimset('TolX',1e-8);
dH_ram = 0.4;
dH_ret = 0.4;
r_th = r_scan;
for ind_r = 1: length(r_th)
    r  = r_th(ind_r);
    tspan(1) = min(-10/r,-1000);
    tspan(2) = max (40/r, 1000);
    fun_ram = @(input) does_it_tip_ram(r,input,H0,tspan,init,opts);
    fun_ret = @(input) does_it_tip_ret(r,input,H0,tspan,init,opts);
    dH_ram = fzero(fun_ram,dH_ram,fzeroopts);
    dH_ret = fzero(fun_ret,dH_ret,fzeroopts);
    dH_ram_th(ind_r) = dH_ram;
    dH_ret_th(ind_r) = dH_ret;
    disp(ind_r)
end
%% BI
dH_ram = 0.4;
r = inf;
tspan(1) = min(-10/r,-1000);
tspan(2) = max (40/r, 1000);
fun_BI = @(input) does_it_tip_ram(r,input,H0,tspan,init,opts);
dH_BI = fzero(fun_BI,dH_ram);

%% testing 
dH = 0.4055556;
r = 1.456^-6;

H  = @(tt)((H0 + dH * sech(r*tt))*(tt<0) + dH*(tt>=0));
% H  = @(tt)(H0 + dH * sech(r*tt));
odefun = @(tt,x)BoxModel_2DH_IVP(tt,x,H);
[time,var_test]   = ode45(odefun,[-10000000, 1000000],init,opts);
% AMOC_off = [lo_bran1(dH), lo_bran2(dH)]
% norm(var_test(end,:)- AMOC_off) < EPS
figure(2)
plot(time,var_test(:,2))
%% plotting 
figure(1);
cla
hold on
plott = pcolor(dH_scan,r_scan,tipping_mat);
plott.LineStyle = 'none';
plott.FaceAlpha = 0.3;

plot(dH_ram_th,r_scan,'--k','LineWidth',1.5)
plot(dH_ret_th,r_scan,'-k','LineWidth',1.5)

colormap([1 0 0; 0 1 0; 1 1 1]);
set(gca, 'YScale', 'log')
set(gcf, 'renderer', 'painters')
axis([0.35 0.45 0 0.1])

%% fuction

% tipping threshold functions 
% ramp shift 
function output = does_it_tip_ram(r,dH,H0,tspan,init, opts)
% output -1 no tipping, or 1 tipping
H  = @(tt)((H0 + dH * sech(r*tt))*(tt<0) + dH*(tt>=0));
odefun = @(tt,x)BoxModel_2DH_IVP(tt,x,H);
[~,var]   = ode45(odefun,tspan,init,opts);
if ( var(end,1) < -0.1)
    output = -1;
else
    output = 1;
end
end

% return shift
function output = does_it_tip_ret(r,dH,H0,tspan,init, opts)
% output -1 no tipping, or 1 tipping
H  = @(tt)(H0 + dH * sech(r*tt));
odefun = @(tt,x)BoxModel_2DH_IVP(tt,x,H);
[~,var]   = ode45(odefun,tspan,init,opts);
if (var(end,1) < -0.1)
    output = -1;
else
    output = 1;
end
end