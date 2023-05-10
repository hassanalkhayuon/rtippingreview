% phase digrams of the AMOC 3 box model showing R-tipping with monotonic fresh 
% water hosing 
warning off
clc
clear
set(0,'defaulttextInterpreter','latex');


opts = odeset('RelTol',1e-5,'AbsTol',1e-10);
load('data/1parBD_AMOC.mat')

up_bran1  = @(input)interp1([us(1:end-1,1);uu(:,1)],...
    [us(1:end-1,2);uu(:,2)],input,'spline');
lo_bran1 = @(input)interp1(ds(:,1),ds(:,2),...
     input,'spline');
un_bran1 = @(input)interp1(du(:,1),du(:,2),...
      input,'spline');

up_bran2  = @(input)interp1([us(1:end-1,1);uu(:,1)],...
    [us(1:end-1,3);uu(:,3)],input);
lo_bran2 = @(input)interp1(ds(:,1),ds(:,3),...
     input);
un_bran2 = @(input)interp1(du(:,1),du(:,3),...
      input);

H0 = 0;
dH = .365;

r_slow    = 0.005;
r_crit1   = 0.0164475257911544780335;
r_fast    = 0.017;
r_crit2   = 0.0178041868523256177;
r_faster  = 1;


%par = [H0,r_slow,dH];
% H_slow =@(tt) ( dH*( 1+tanh(r_slow*tt)) /2 ) ;
% H_crit1 =@(tt) ( dH*( 1+tanh(r_crit1*tt)) /2 ) ;
% H_fast =@(tt) ( dH*( 1+tanh(r_fast*tt)) /2 ) ;
% H_crit2 =@(tt) ( dH*( 1+tanh(r_crit2*tt)) /2 ) ;
% H_faster =@(tt) ( dH*( 1+tanh(r_faster*tt)) /2 ) ;


H_slow   =@(tt) ((H0 + dH * sech(r_slow*tt))*(tt<0) + (H0 + dH)*(tt>=0));
H_crit1  =@(tt) ((H0 + dH * sech(r_crit1*tt))*(tt<0) + (H0 + dH)*(tt>=0));
H_fast   =@(tt) ((H0 + dH * sech(r_fast*tt))*(tt<0) + (H0 + dH)*(tt>=0)); 
H_crit2  =@(tt) ((H0 + dH * sech(r_crit2*tt))*(tt<0) + (H0 + dH)*(tt>=0));
H_faster =@(tt) ((H0 + dH * sech(r_faster*tt))*(tt<0) + (H0 + dH)*(tt>=0));

tspan_slow        = [-2000, 4000];
tspan_crit1       = [-2000, 4000];
tspan_fast        = [-3000, 4000];
tspan_crit2       = [-3000, 4000];
tspan_faster      = [-1200, 4000];

init         = [up_bran1(H0),up_bran2(H0)]; % define your initial point

% slow trajectory 
odefun_slow       = @(tt,x)BoxModel_2DH_IVP(tt,x,H_slow);
[time_slow,x_slow]        = ode45(odefun_slow,tspan_slow ,init,opts);

% tipping threshold trajectory 1
odefun_crit1       = @(tt,x)BoxModel_2DH_IVP(tt,x,H_crit1);
[time_crit1,x_crit1]        = ode45(odefun_crit1,tspan_crit1 ,init,opts);

% fast trajectory
odefun_fast       = @(tt,x)BoxModel_2DH_IVP(tt,x,H_fast);
[time_fast,x_fast]        = ode45(odefun_fast,tspan_fast ,init,opts);

% tipping threshold trajectory 1
odefun_crit2       = @(tt,x)BoxModel_2DH_IVP(tt,x,H_crit2);
[time_crit2,x_crit2]        = ode45(odefun_crit2,tspan_crit2,init,opts);

% faster trajectory
odefun_faster     = @(tt,x)BoxModel_2DH_IVP(tt,x,H_faster);
[time_faster,x_faster]        = ode45(odefun_faster,tspan_faster ,init,opts);

%% Coloring the other basin of attraction 
% res = 500;
% x_scan      = linspace(-0.23,0.1,res);
% y_scan      = linspace(0,0.35,res);
% HH          = @(tt)dH;
% ON_EQ       = [up_bran1(HH(0)),up_bran2(HH(0))];
% tspan_basin = [0 3000];
% colormat    = ones(res,res);
% 
% for ind_x = 1:res
%     parfor ind_y = 1:res
%         init_basin         = [x_scan(ind_x);y_scan(ind_y)];
%         odefun_basin       = @(tt,x)BoxModel_2DH_IVP(tt,x,HH);
%         [~,x_basin]        = ode45(odefun_basin,tspan_basin,init_basin);
%         if norm(x_basin(end,:) - ON_EQ) < 0.01
%             colormat(ind_y,ind_x) = 0;
%         end
%     end
%     disp(ind_x);
% end
load('data/basin_Hmax.mat')
% grayscal = 1:-0.01:0;
% map = [grayscal',grayscal', grayscal'];
% PCOLOR = pcolor(x_scan,y_scan,colormat);
% PCOLOR.FaceAlpha = 0.1;
% PCOLOR.LineStyle = 'none';
% colormap(map)

%%

% upplim = 0.0323746;
% lowlim = -0.0824422;
% 
% threshold = @(xx)interp1(x_crit(:,1),x_crit(:,2),xx);
% 
% 
% 
% resl = 1000;
% 
% x_scan = linspace(-0.3,0.05,resl);
% y_scan = linspace(0,0.3,resl);
% colormat = zeros(resl,resl);
% 
% 
% for ind_x = 1:resl
%     
%     X  = x_scan(ind_x);
%     
%     % Region 1
%     
%     if X < lowlim 
%         parfor ind_y = 1:resl
%                 colormat(ind_y,ind_x) = 0;
%         end
%     end
%         
%     %     Region 2
% 
%     if X > upplim
%         parfor ind_y = 1:resl
%             Y  = y_scan(ind_y);
%             if Y>=0.14345
%                 colormat(ind_y,ind_x) = 1;
%             end
%         end
%     end
%     
%     %     Region 3
% 
%     if (X <= upplim && X >= lowlim)
%         Ther = threshold(X);
%         parfor ind_y = 1:resl
%             Y  = y_scan(ind_y);
%             if Y>=Ther
%                 colormat(ind_y,ind_x) = 1;
%             end
%         end
%     end
%     
%     disp(ind_x);
% end
%% trajectories MIN

% HH          = @(tt)H0;
% tspan_basin = [0 350];
% 
% % The ON inits 
% init_basin        = [-0.003;0.06];
% % init_basin        = [-0.056;0.191];
% % init_basin        = [-0.069;0.3];
% % init_basin        = [0.07;0.03];
% 
% %The off inits
% 
% %init_basin        = [-0.217;0.022]; Extra
% %init_basin        = [-0.138;0.0182];
% 
% 
% % init_basin        = [-0.21,0.33];
% % init_basin        = [-0.11,0.2];
% % init_basin        = [-0.06,0.033];
% % init_basin        = [-0.107,0.03086];
% % init_basin        = [-0.185,0.029];
% 
% odefun_basin       = @(tt,x)BoxModel_2DH_IVP(tt,x,HH);
% [~,x_var]        = ode45(odefun_basin,tspan_basin,init_basin,opts);
% 
% %plot(init_basin(1),init_basin(2),'.r','MarkerSize',20)
% plot(x_var(:,1),x_var(:,2),'-r','LineWidth',2)

%% trajectories MAX

% HH          = @(tt)dH;
% tspan_basin = [0 2000];
% 
% % The ON inits 
%  init_basin        = [-0.065;0.246];
% 
%  % The Off inits 
% init_basin        = [0.037,0.3];
% init_basin        = [-0.01507,0.0231];
% init_basin        = [-0.1,0.281];
% 
% 
% odefun_basin       = @(tt,x)BoxModel_2DH_IVP(tt,x,HH);
% [~,x_var]        = ode45(odefun_basin,tspan_basin,init_basin,opts);
% 
% plot(init_basin(1),init_basin(2),'.r','MarkerSize',20)
% %plot(x_var(:,1),x_var(:,2),'-r','LineWidth',2)

%% Plot 

% new colors 

red = [0.86,0.02,0.05];
yellow = [221,170,51]/255;
green = [0.31,0.70,0.40];
blue = [0.10,0.40,0.69];

H_branch = linspace(H0,dH,100);

On_branch = [up_bran1(H_branch); up_bran2(H_branch)];


e_on_max = [up_bran1(dH); up_bran2(dH)];
e_on_min = [up_bran1(H0); up_bran2(H0)];

e_off_max = [lo_bran1(dH); lo_bran2(dH)];
e_off_min = [lo_bran1(H0); lo_bran2(H0)];

e_s_max = [un_bran1(dH); un_bran2(dH)];
e_s_min = [un_bran1(H0); un_bran2(H0)];

[~,w_s_max1] = ode45 (@(t,x)BoxModel_2DH_IVP(t,x,@(t)dH), ...
    [10000 0],e_s_max+[0; 0.001]);
[~,w_s_max2] = ode45 (@(t,x)BoxModel_2DH_IVP(t,x,@(t)dH), ...
    [5000 0],e_s_max-[0; 0.001]);

[~,w_s_min1] = ode45 (@(t,x)BoxModel_2DH_IVP(t,x,@(t)H0), ...
    [10000 0],e_s_min+[0; 0.001]);
[~,w_s_min2]= ode45 (@(t,x)BoxModel_2DH_IVP(t,x,@(t)H0), ...
    [10000 0],e_s_min-[0; 0.001]);

figure;
cla
box on
disableDefaultInteractivity(gca)
hold on

grayscal = 1:-0.01:0;
map = [grayscal',grayscal', grayscal'];
PCOLOR = pcolor(x_scan,y_scan,colormat);
PCOLOR.FaceAlpha = 0.1;
PCOLOR.LineStyle = 'none';
colormap(map)


plot(...
    e_s_max(1),e_s_max(2),'+',...
    'MarkerSize',10,'Color','k','LineWidth',2)
plot(...
    w_s_min1(:,1),w_s_min1(:,2),'--',...
    w_s_min2(:,1),w_s_min2(:,2),'--',...
    'Color',[0.7 0.7,0.7],'LineWidth',1)
figure; 
cla
hold on
plot(...
    w_s_max1(:,1),w_s_max1(:,2),'--',...
    w_s_max2(:,1),w_s_max2(:,2),'--',...
    'Color',[0.7 0.7,0.7],'LineWidth',1)


% critical trajectory 1
 plot(x_crit1(1:300,1),x_crit1(1:300,2),'-',...
     'Color',blue,'LineWidth',2)

% slow trajectory 
plot(x_slow(:,1),x_slow(:,2),'-',...
    'Color',green,'LineWidth',2)
% 
% % % fast trajectory
% plot(x_fast(:,1),x_fast(:,2),'-',...
%     'Color',[.4 .2 .4],'LineWidth',2)
% % 
% % % critical trajectory 2
% % 
% plot(x_crit2(:,1),x_crit2(:,2),'-',...
%      'Color',[1 .4 .15],'LineWidth',2)
% 
% even faster trajectories
plot(x_faster(:,1),x_faster(:,2),'-',...
    'Color',red,'LineWidth',2)


plot(...
    e_on_min(1),e_on_min(2),'.',...
    'MarkerSize',30,'Color','k')

plot(...
    e_off_min(1),e_off_min(2),'.',...
    'MarkerSize',30,'Color','k')

plot(...
    e_on_max(1),e_on_max(2),'.',...
    e_off_max(1),e_off_max(2),'.', ...
    'MarkerSize',30,'Color','k')

plot(...
    On_branch(1,:),On_branch(2,:),':',...
    'Color','k','LineWidth',1.5);

plot(...
    e_s_min(1),e_s_min(2),'+',...
    'MarkerSize',10,'Color',[0.7,0.7,0.7],'LineWidth',2)

box on
axis([-0.23 0.1 0 0.35])
set(gca,'FontSize',16)

%% plotting time series

Ts = 7.919; %C
T0 = 3.87; %C
lambda = 1.62e7; %m^6kg^-1s^-1
alpha = 0.12; %kg m^-3 C^-1
beta = 790.0;%kg m^-3 
mu = 22e-8; %deg s m^-3

SN_crit =  x_crit1(:,1);
SN_slow =  x_slow(:,1);
SN_fast =  x_faster(:,1);

S0 = 0.035; % dimensionless

SS = (0.034427-S0).*100;

Q_crit = lambda.*(alpha.*(Ts-T0)+beta.*(SN_crit./100-SS./100))/(1+lambda*alpha*mu)*10^-6;
Q_fast = lambda.*(alpha.*(Ts-T0)+beta.*(SN_fast./100-SS./100))/(1+lambda*alpha*mu)*10^-6;
Q_slow = lambda.*(alpha.*(Ts-T0)+beta.*(SN_slow./100-SS./100))/(1+lambda*alpha*mu)*10^-6;

% Q units is Sv (1 million cubic metres per second)

figure(2)
cla
hold on

% critical rate
plot(time_crit1,Q_crit,'-',...
    'Color',yellow,'LineWidth',2)

% fast 

%even faster trajectories
plot(time_faster,Q_fast,'-',...
    'Color',red,'LineWidth',2)

% slow 

plot(time_slow,Q_slow,'-',...
    'Color',blue,'LineWidth',2)

xlim([-1000, 3090])

% figure
% hold on
fplot(H_slow, [-1000, 4000],'Color',blue,'LineWidth',2)
fplot(H_crit1, [-1000, 4000],'Color',yellow,'LineWidth',2)
fplot(H_faster, [-1000, 4000],'Color',red,'LineWidth',2)
% ylim([0 0.4])




