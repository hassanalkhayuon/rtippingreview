% A movie of the phase digrams of the AMOC 3 box model showing R-tipping
% with monotonic fresh water hosing

warning off
clc
clear
set(0,'defaulttextInterpreter','latex');

opts = odeset('RelTol',1e-5,'AbsTol',1e-10);
load('data/1parBD_AMOC.mat')
ds(:,3) = ds(:,3) - 0.0029;
%% equilibria branches and parametrs

H0 = 0;
dH = .365;

Ts = 7.919; %C
T0 = 3.87; %C
lambda = 1.62e7; %m^6kg^-1s^-1
alpha = 0.12; %kg m^-3 C^-1
beta = 790.0;%kg m^-3
mu = 22e-8; %deg s m^-3

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
%% Three trajectories with dynamic hosing

r_superslow  = 0.0001;
r_slow  = 0.005;
r_crit  = 0.0164491186245;
r_fast  = 0.05;

H_superslow  = @(tt) ((H0 + dH * sech(r_superslow*tt))*(tt<0) + (H0 + dH)*(tt>=0));
H_slow  = @(tt) ((H0 + dH * sech(r_slow*tt))*(tt<0) + (H0 + dH)*(tt>=0));
H_crit  = @(tt) ((H0 + dH * sech(r_crit*tt))*(tt<0) + (H0 + dH)*(tt>=0));
H_fast  = @(tt) ((H0 + dH * sech(r_fast*tt))*(tt<0) + (H0 + dH)*(tt>=0));

tspan_slow   = -2000:5:3200;
tspan_crit   = -2000:5:3200;
tspan_fast   = -2000:5:3200;

init = [up_bran1(H0),up_bran2(H0)]; % initial point (the on-state at H0)

% slow trajectory
odefun_slow   = @(tt,x)BoxModel_2DH_IVP(tt,x,H_slow);
[~,x_slow]    = ode45(odefun_slow,tspan_slow ,init,opts);

% tipping threshold trajectory 1
odefun_crit   = @(tt,x)BoxModel_2DH_IVP(tt,x,H_crit);
[~,x_crit]    = ode45(odefun_crit,tspan_crit ,init,opts);


% fast trajectory
odefun_fast      = @(tt,x)BoxModel_2DH_IVP(tt,x,H_fast);
[time,x_fast]    = ode45(odefun_fast,tspan_fast ,init,opts);

% to find AMOC strength

SN_crit =  x_crit(:,1);
SN_slow =  x_slow(:,1);
SN_fast =  x_fast(:,1);

S0 = 0.035; % dimensionless

SS = (0.034427-S0).*100;

Q_crit = lambda.*(alpha.*(Ts-T0)+beta.*(SN_crit./100-SS./100))/(1+lambda*alpha*mu)*10^-6;
Q_fast = lambda.*(alpha.*(Ts-T0)+beta.*(SN_fast./100-SS./100))/(1+lambda*alpha*mu)*10^-6;
Q_slow = lambda.*(alpha.*(Ts-T0)+beta.*(SN_slow./100-SS./100))/(1+lambda*alpha*mu)*10^-6;

%% Plots and movie

% new colors
red    = [0.86,0.02,0.05];
yellow = [0.87 0.67 0.20];
green  = [0.31,0.70,0.40];
blue   = [0.10,0.40,0.69];
fainted_red   = [0.86,0.02,0.05,0.4];
fainted_green = [0.31,0.70,0.40,0.4];
fainted_blue  = [0.10,0.40,0.69,0.4];

for ind = 1:length(time)
    H_C_slow (ind)  = H_slow(time(ind));
    SN_on_QS_slow (ind)  = up_bran1(H_slow(time(ind)));
    SN_off_QS_slow (ind) = lo_bran1(H_slow(time(ind)));

    H_C_crit (ind)  = H_crit(time(ind));
    SN_on_QS_crit (ind)  = up_bran1(H_crit(time(ind)));
    SN_off_QS_crit (ind) = lo_bran1(H_crit(time(ind)));
    
    H_C_fast (ind)  = H_fast(time(ind));
    SN_on_QS_fast (ind)  = up_bran1(H_fast(time(ind)));
    SN_off_QS_fast (ind) = lo_bran1(H_fast(time(ind)));
    

end

Q_on_slow  = lambda.*(alpha.*(Ts-T0)+beta.*(SN_on_QS_slow./100-SS./100))/(1+lambda*alpha*mu)*10^-6;
Q_off_slow = lambda.*(alpha.*(Ts-T0)+beta.*(SN_off_QS_slow./100-SS./100))/(1+lambda*alpha*mu)*10^-6;

Q_on_fast  = lambda.*(alpha.*(Ts-T0)+beta.*(SN_on_QS_fast./100-SS./100))/(1+lambda*alpha*mu)*10^-6;
Q_off_fast = lambda.*(alpha.*(Ts-T0)+beta.*(SN_off_QS_fast./100-SS./100))/(1+lambda*alpha*mu)*10^-6;

Q_on_crit  = lambda.*(alpha.*(Ts-T0)+beta.*(SN_on_QS_crit./100-SS./100))/(1+lambda*alpha*mu)*10^-6;
Q_off_crit = lambda.*(alpha.*(Ts-T0)+beta.*(SN_off_QS_crit./100-SS./100))/(1+lambda*alpha*mu)*10^-6;

%% Slow forcing movie 

movie_name_s = 'AMOC_Rtipping_temp';
V_slow = VideoWriter(movie_name_s,'MPEG-4');
V_slow.FrameRate = 30;
V_slow.Quality = 100;

open(V_slow);
end_of_run = find(time >= 2000,1);
for ind = 200:end_of_run

    % the on- off- and saddle- states
    S_on  = [up_bran1(H_slow(time(ind)));...
        up_bran2(H_slow(time(ind)))];
    S_off = [lo_bran1(H_slow(time(ind)));...
        lo_bran2(H_slow(time(ind)))];
    S_sad =  [un_bran1(H_slow(time(ind)));...
        un_bran2(H_slow(time(ind)))];

    % tipping boundary -- the stable manifold on S_sad
    HH = @(tt)H_slow(time(ind));
    [~,w_s] = ode45 (@(t,x)BoxModel_2DH_IVP(t,x,HH), ...
        7000:-10:0,S_sad-[0; 0.001]);

    figure(10);
    clf

    % first subplot (forcing time series)
    h1 = subplot(7,13,[1:6,14:19,27:32]);
    hold on

    % slow forcing 
    plot(...
        time,H_C_slow,...
        'LineWidth',2,'Color',fainted_green)
    plot(...
        time(1:ind), H_C_slow(1:ind),...
        'LineWidth',2,'Color',green)
    plot(...
        time(ind),H_C_slow(ind),'.',...
        'MarkerSize',25,'Color',green)   
    box off
    xlim([-1000 3200])
    ylim([0 0.4])
    xticks([-1000 0 1000 2000 3000])
    xticklabels({})
    yticks([0 0.2 0.4])
    yticklabels([0.0 0.2 0.4])
    ylabel(['Freshwater ';'Hosing [Sv]']);

    % second subplot (salinity time series)
    h2 = subplot(7,13,[53:58,66:71,79:84]);
    hold on

    % QSE
    plot(...
        time,Q_on_slow,'--k',...
        time,Q_off_slow,'--k',...
        'LineWidth',1,'Color',[0.7 0.7 0.7])

    %time series 
    plot(...
        time(1:ind),Q_slow(1:ind),'Color',green,...
        'LineWidth',2)
    plot(...
        time(ind),Q_slow(ind),'.','Color',green,...
        'MarkerSize',25)

    xlim([-1000 3200])
    ylim([-10 15])
    xticks([-1000 0 1000 2000 3000])
    xticklabels([0 1000 2000 3000 4000])
    yticks([-10 0 10])
    yticklabels([-10 0 10])
    
    xlabel('time [years]')
    ylabel(['AMOC         ';'Strength [Sv]']);

    % third subplot (phase space)
    h3 = subplot(7,13,[8:13,21:26,34:39,60:65,73:78,86:91]);
    hold on
    
    % slow forcing
    if time(ind)<=-25
    plot(...
        w_s(:,1),w_s(:,2),'-',...
        'LineWidth',2,'Color',0.75*[1 1 1])
    else
        plot(...
            w_s(500:end,1),w_s(500:end,2),'-',...
            'LineWidth',2,'Color',0.75*[1 1 1])
    end
    plot(...
        x_slow(1:ind,1),x_slow(1:ind,2),...
        'LineWidth',2,'Color',green)
    plot(...
        x_slow(1,1),x_slow(1,2),'k.',...
        'MarkerSize',30)
    plot(...
        x_slow(ind,1),x_slow(ind,2),'.',...
        'MarkerSize',25,'Color',green)

    plot(...
        S_on(1),S_on(2),'k.',...
        'MarkerSize',30)
    plot(...
        S_off(1),S_off(2),'k.',...
        'MarkerSize',30)
    plot(...
        S_sad(1),S_sad(2),'ko',...
        'MarkerSize',6.5)

    xlabel('North Atlantic Salinity [psu]')
    ylabel('Tropical Atlantic Salinity [psu]')
    
    set(h1,'FontSize',20)
    set(h1,'Position',[0.09 0.6 0.43  0.33])

    set(h2,'FontSize',20)
    set(h2,'Position',[0.09 0.12 0.43 0.33])

    set(h3,'FontSize',20)
    set(h3,'Position',[0.62 0.12 0.35  0.82])


    xlim([-0.25 0.05]);
    ylim([0 0.3]);
    xticks([-0.2 -0.1 0])
    xticklabels([33 34 35])
    yticks([0 0.1 0.2 0.3])
    yticklabels([35 36 37 38])
    box on

    left   =  .3e3;
    bottom =  0.09e3;
    width  =  1.6e3;
    height =  .7e3;

    set(gcf,'position', [left bottom width height])
    set(gcf,'color','w');  
    set(gcf, 'renderer', 'painters')
    writeVideo(V_slow,getframe(gcf));
%     drawnow
end
close(V_slow);

%% fast forcing movie + slow focing in the background

movie_name_f = 'AMOC_Rtipping_f';
V_fast = VideoWriter(movie_name_f,'MPEG-4');
V_fast.FrameRate = 30;
V_fast.Quality = 100;

open(V_fast);
end_of_run = find(time >= 2000,1);
for ind = 200:end_of_run

    % the on- off- and saddle- states
    S_on  = [up_bran1(H_fast(time(ind)));...
        up_bran2(H_fast(time(ind)))];
    S_off = [lo_bran1(H_fast(time(ind)));...
        lo_bran2(H_fast(time(ind)))];
    S_sad =  [un_bran1(H_fast(time(ind)));...
        un_bran2(H_fast(time(ind)))];

    % tipping boundary -- the stable manifold on S_sad
    HH = @(tt)H_fast(time(ind));
    [~,w_f] = ode45 (@(t,x)BoxModel_2DH_IVP(t,x,HH), ...
        7000:-10:0,S_sad-[0; 0.001]);

    figure(10);
    clf

    % first subplot (forcing time series)
    h1 = subplot(7,13,[1:6,14:19,27:32]);
    hold on
    % slow forcing 
    plot(...
        time,H_C_slow,...
        'LineWidth',2,'Color',fainted_green)
    % fast forcing
    plot(...
        time,H_C_fast,...
        'LineWidth',2,'Color',fainted_red)
    plot(...
        time(1:ind), H_C_fast(1:ind),...
        'LineWidth',2,'Color',red)
    plot(...
        time(ind),H_C_fast(ind),'.',...
        'MarkerSize',25,'Color',red)   
    box off
    xlim([-1000 3200])
    ylim([0 0.4])
    xticks([-1000 0 1000 2000 3000])
    xticklabels({})
    yticks([0 0.2 0.4])
    yticklabels([0.0 0.2 0.4])
    ylabel(['Freshwater ';'Hosing [Sv]']);

    % second subplot (salinity time series)
    h2 = subplot(7,13,[53:58,66:71,79:84]);
    hold on

    % QSE
    plot(...
        time,Q_on_fast,'--k',...
        time,Q_off_fast,'--k',...
        'LineWidth',1,'Color',[0.75 0.75 0.75])

    %time series 
    plot(...
        time,Q_slow,'Color',fainted_green,...
        'LineWidth',2)
    plot(...
        time(end),Q_slow(end),'.','Color', ...
        fainted_green,'MarkerSize',25)
    plot(...
        time(1:ind),Q_fast(1:ind),'Color',red,...
        'LineWidth',2)
    plot(...
        time(ind),Q_fast(ind),'.','Color',red,...
        'MarkerSize',25)

    xlim([-1000 3200])
    ylim([-10 15])
    xticks([-1000 0 1000 2000 3000])
    xticklabels([0 1000 2000 3000 4000])
    yticks([-10 0 10])
    yticklabels([-10 0 10])
    
    xlabel('time [years]')
    ylabel(['AMOC         ';'Strength [Sv]']);

    % third subplot (phase space)
    h3 = subplot(7,13,[8:13,21:26,34:39,60:65,73:78,86:91]);
    hold on
    
    % fast forcing
    if time(ind)<=0
        plot(...
            w_f(:,1),w_f(:,2),'-',...
            'LineWidth',2,'Color',0.75*[1 1 1])
    else
        plot(...
            w_f(500:end,1),w_f(500:end,2),'-',...
            'LineWidth',2,'Color',0.75*[1 1 1])
    end
    plot(...
        x_slow(:,1),x_slow(:,2),...
        'LineWidth',2,'Color',fainted_green)
    plot(...
        x_slow(end,1),x_slow(end,2),'.',...
        'Color',green,'MarkerSize',25)
    plot(...
        x_fast(1:ind,1),x_fast(1:ind,2),...
        'LineWidth',2,'Color',red)
    plot(...
        x_fast(1,1),x_fast(1,2),'.k',...
        'MarkerSize',30)
    plot(...
        x_fast(ind,1),x_fast(ind,2),'.',...
        'MarkerSize',25,'Color',red)

    plot(...
        S_on(1),S_on(2),'k.',...
        'MarkerSize',30)
    plot(...
        S_off(1),S_off(2),'k.',...
        'MarkerSize',30)
    plot(...
        S_sad(1),S_sad(2),'ko',... 
        'MarkerSize',6.5)

    xlabel('North Atlantic Salinity [psu]')
    ylabel('Tropical Atlantic Salinity [psu]')
    
    set(h1,'FontSize',20)
    set(h1,'Position',[0.09 0.6 0.43  0.33])

    set(h2,'FontSize',20)
    set(h2,'Position',[0.09 0.12 0.43 0.33])

    set(h3,'FontSize',20)
    set(h3,'Position',[0.62 0.12 0.35  0.82])


    xlim([-0.25 0.05]);
    ylim([0 0.3]);
    xticks([-0.2 -0.1 0])
    xticklabels([33 34 35])
    yticks([0 0.1 0.2 0.3])
    yticklabels([35 36 37 38])
    box on

    left   =  .3e3;
    bottom =  0.09e3;
    width  =  1.6e3;
    height =  .7e3;

    set(gcf,'position', [left bottom width height])
    set(gcf,'color','w');  
    set(gcf, 'renderer', 'painters')
    writeVideo(V_fast,getframe(gcf));
end
close(V_fast)

%% critical forcing movie 

movie_name_c = 'AMOC_Rtipping_c';
V_crit = VideoWriter(movie_name_c,'MPEG-4');
V_crit.FrameRate = 30;
V_crit.Quality = 100;

open(V_crit);
for ind = 200:length(time)

    % the on- off- and saddle- states
    S_on  = [up_bran1(H_crit(time(ind)));...
        up_bran2(H_crit(time(ind)))];
    S_off = [lo_bran1(H_crit(time(ind)));...
        lo_bran2(H_crit(time(ind)))];
    S_sad =  [un_bran1(H_crit(time(ind)));...
        un_bran2(H_crit(time(ind)))];

    % tipping boundary -- the stable manifold on S_sad
    HH = @(tt)H_crit(time(ind));
    [~,w_c] = ode45 (@(t,x)BoxModel_2DH_IVP(t,x,HH), ...
        7000:-10:0,S_sad-[0; 0.001]);

    figure(10);
    clf

    % first subplot (forcing time series)
    h1 = subplot(7,13,[1:6,14:19,27:32]);
    hold on
    % slow forcing 
    plot(...
        time,H_C_slow,...
        'LineWidth',2,'Color',fainted_green)
    % fast forcing
    plot(...
        time,H_C_fast,...
        'LineWidth',2,'Color',fainted_red)  
    % critical forcing 
    plot(...
        time,H_C_crit,...
        'LineWidth',2,'Color',fainted_blue)
    plot(...
        time(1:ind), H_C_crit(1:ind),...
        'LineWidth',2,'Color',blue)
    plot(...
        time(ind),H_C_crit(ind),'.',...
        'MarkerSize',25,'Color',blue)

    box off
    xlim([-1000 3200])
    ylim([0 0.4])
    xticks([-1000 0 1000 2000 3000])
    xticklabels({})
    yticks([0 0.2 0.4])
    yticklabels([0.0 0.2 0.4])
    ylabel(['Freshwater ';'Hosing [Sv]']);

    % second subplot (salinity time series)
    h2 = subplot(7,13,[53:58,66:71,79:84]);
    hold on

    % QSE
    plot(...
        time,Q_on_crit,'--k',...
        time,Q_off_crit,'--k',...
        'LineWidth',1,'Color',[0.7 0.7 0.7])

    %time series 
    plot(...
        time,Q_slow,'Color',fainted_green,...
        'LineWidth',2)
    plot(...
        time(end),Q_slow(end),'.','Color', ...
        fainted_green,'MarkerSize',25)
    plot(...
        time,Q_fast,'Color',fainted_red,...
        'LineWidth',2)
    plot(...
        time(end),Q_fast(end),'.','Color', ...
        fainted_red,'MarkerSize',25)
    plot(...
        time(1:ind),Q_crit(1:ind),'Color',blue,...
        'LineWidth',2)
    plot(...
        time(ind),Q_crit(ind),'.','Color',blue,...
        'MarkerSize',25)

    xlim([-1000 3200])
    ylim([-10 15])
    xticks([-1000 0 1000 2000 3000])
    xticklabels([0 1000 2000 3000 4000])
    yticks([-10 0 10])
    yticklabels([-10 0 10])
    
    xlabel('time [years]')
    ylabel(['AMOC         ';'Strength [Sv]']);

    % third subplot (phase space)
    h3 = subplot(7,13,[8:13,21:26,34:39,60:65,73:78,86:91]);
    hold on
    
    % critical forcing

    if time(ind)<=-5
        plot(...
            w_c(:,1),w_c(:,2),'-','Color',0.75.*[1 1 1],...
            'LineWidth',2)
    else
        plot(...
            w_c(500:end,1),w_c(500:end,2),'-','Color',0.75.*[1 1 1],...
            'LineWidth',2)
    end

    plot(...
        x_slow(:,1),x_slow(:,2),'Color',fainted_green,...
        'LineWidth',2)
    plot(...
        x_slow(end,1),x_slow(end,2),'.','Color',fainted_green,...
        'MarkerSize',25)
    plot(...
        x_fast(end,1),x_fast(end,2),'.','Color',fainted_red,...
        'MarkerSize',25)
    plot(x_fast(:,1),x_fast(:,2),'Color',fainted_red,...
        'LineWidth',2)    
    plot(...
        x_crit(1:ind,1),x_crit(1:ind,2),'-',...
        'Color',blue,'LineWidth',2)
    plot(x_crit(1,1),x_crit(1,2),'.k',...
        'MarkerSize',30)
    plot(...
        x_crit(ind,1),x_crit(ind,2),'.',...
        'Color',blue,'MarkerSize',25)


    plot(...
        S_on(1),S_on(2),'k.',...
        'MarkerSize',30)
    plot(...
        S_off(1),S_off(2),'k.',...
        'MarkerSize',30)
    plot(...
        S_sad(1),S_sad(2),'ko',... 
        'MarkerSize',6.5)

    xlabel('North Atlantic Salinity [psu]')
    ylabel('Tropical Atlantic Salinity [psu]')
    
    set(h1,'FontSize',20)
    set(h1,'Position',[0.09 0.6 0.43  0.33])

    set(h2,'FontSize',20)
    set(h2,'Position',[0.09 0.12 0.43 0.33])

    set(h3,'FontSize',20)
    set(h3,'Position',[0.62 0.12 0.35  0.82])

    xlim([-0.25 0.05]);
    ylim([0 0.3]);
    xticks([-0.2 -0.1 0])
    xticklabels([33 34 35])
    yticks([0 0.1 0.2 0.3])
    yticklabels([35 36 37 38])
    box on

    left   =  .3e3;
    bottom =  0.09e3;
    width  =  1.6e3;
    height =  .7e3;

    set(gcf,'position', [left bottom width height])
    set(gcf,'color','w');  
    set(gcf, 'renderer', 'painters')
    writeVideo(V_crit,getframe(gcf));
%     drawnow
end
close(V_crit)
%% the full movie 
movie_name_full = 'AMOC_Rtipping_full';
V_full = VideoWriter(movie_name_full,'MPEG-4');
V_full.FrameRate = 30;
V_full.Quality = 100;

open(V_full);

% first part (slow only)
for ind = 200:length(time)
% the on- off- and saddle- states
    S_on  = [up_bran1(H_slow(time(ind)));...
        up_bran2(H_slow(time(ind)))];
    S_off = [lo_bran1(H_slow(time(ind)));...
        lo_bran2(H_slow(time(ind)))];
    S_sad =  [un_bran1(H_slow(time(ind)));...
        un_bran2(H_slow(time(ind)))];

    % tipping boundary -- the stable manifold on S_sad
    HH = @(tt)H_slow(time(ind));
    [~,w_s] = ode45 (@(t,x)BoxModel_2DH_IVP(t,x,HH), ...
        7000:-10:0,S_sad-[0; 0.001]);

    figure(10);
    clf
    annotation('textbox',...
    [0.412 0.94 0.266 0.06],...
    'String',{'Slow freshwater hosing'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',20);


    % first subplot (forcing time series)
    h1 = subplot(7,13,[1:6,14:19,27:32]);
    hold on

    % slow forcing 
    plot(...
        time,H_C_slow,...
        'LineWidth',2,'Color',fainted_green)
    plot(...
        time(1:ind), H_C_slow(1:ind),...
        'LineWidth',2,'Color',green)
    plot(...
        time(ind),H_C_slow(ind),'.',...
        'MarkerSize',25,'Color',green)   
    box off
    xlim([-1000 3200])
    ylim([0 0.4])
    xticks([-1000 0 1000 2000 3000])
    xticklabels({})
    yticks([0 0.2 0.4])
    yticklabels([0.0 0.2 0.4])
    ylabel(['Freshwater ';'Hosing [Sv]']);

    % second subplot (salinity time series)
    h2 = subplot(7,13,[53:58,66:71,79:84]);
    hold on

    % QSE
    plot(...
        time,Q_on_slow,'--k',...
        time,Q_off_slow,'--k',...
        'LineWidth',1,'Color',[0.7 0.7 0.7])

    %time series 
    plot(...
        time(1:ind),Q_slow(1:ind),'Color',green,...
        'LineWidth',2)
    plot(...
        time(ind),Q_slow(ind),'.','Color',green,...
        'MarkerSize',25)

    xlim([-1000 3200])
    ylim([-10 15])
    xticks([-1000 0 1000 2000 3000])
    xticklabels([0 1000 2000 3000 4000])
    yticks([-10 0 10])
    yticklabels([-10 0 10])
    
    xlabel('time [years]')
    ylabel(['AMOC         ';'Strength [Sv]']);

    % third subplot (phase space)
    h3 = subplot(7,13,[8:13,21:26,34:39,60:65,73:78,86:91]);
    hold on
    
    % slow forcing
    if time(ind)<=-25
    plot(...
        w_s(:,1),w_s(:,2),'-',...
        'LineWidth',2,'Color',0.75*[1 1 1])
    else
        plot(...
            w_s(500:end,1),w_s(500:end,2),'-',...
            'LineWidth',2,'Color',0.75*[1 1 1])
    end
    plot(...
        x_slow(1:ind,1),x_slow(1:ind,2),...
        'LineWidth',2,'Color',green)
    plot(...
        x_slow(1,1),x_slow(1,2),'k.',...
        'MarkerSize',30)
    plot(...
        x_slow(ind,1),x_slow(ind,2),'.',...
        'MarkerSize',25,'Color',green)

    plot(...
        S_on(1),S_on(2),'k.',...
        'MarkerSize',30)
    plot(...
        S_off(1),S_off(2),'k.',...
        'MarkerSize',30)
    plot(...
        S_sad(1),S_sad(2),'ko',...
        'MarkerSize',6.5)

    xlabel('North Atlantic Salinity [psu]')
    ylabel('Tropical Atlantic Salinity [psu]')
    
    set(h1,'FontSize',20)
    set(h1,'Position',[0.09 0.59 0.43 0.33])

    set(h2,'FontSize',20)
    set(h2,'Position',[0.09 0.11 0.43 0.33])

    set(h3,'FontSize',20)
    set(h3,'Position',[0.62 0.11 0.35 0.82])


    xlim([-0.25 0.05]);
    ylim([0 0.3]);
    xticks([-0.2 -0.1 0])
    xticklabels([33 34 35])
    yticks([0 0.1 0.2 0.3])
    yticklabels([35 36 37 38])
    box on

    left   =  .3e3;
    bottom =  0.09e3;
    width  =  1.6e3;
    height =  .7e3;

    set(gcf,'position', [left bottom width height])
    set(gcf,'color','w');  
    set(gcf, 'renderer', 'painters')
    writeVideo(V_full,getframe(gcf));
end

% second part (fast and slow)
for ind = 200:length(time)
    % the on- off- and saddle- states
    S_on  = [up_bran1(H_fast(time(ind)));...
        up_bran2(H_fast(time(ind)))];
    S_off = [lo_bran1(H_fast(time(ind)));...
        lo_bran2(H_fast(time(ind)))];
    S_sad =  [un_bran1(H_fast(time(ind)));...
        un_bran2(H_fast(time(ind)))];

    % tipping boundary -- the stable manifold on S_sad
    HH = @(tt)H_fast(time(ind));
    [~,w_f] = ode45 (@(t,x)BoxModel_2DH_IVP(t,x,HH), ...
        7000:-10:0,S_sad-[0; 0.001]);

    figure(10);
    clf
    annotation('textbox',...
        [0.412 0.94 0.266 0.06],...
        'String',{'Fast freshwater hosing'},...
        'LineStyle','none',...
        'Interpreter','latex',...
        'FontSize',20);

    % first subplot (forcing time series)
    h1 = subplot(7,13,[1:6,14:19,27:32]);
    hold on
    % slow forcing 
    plot(...
        time,H_C_slow,...
        'LineWidth',2,'Color',fainted_green)
    % fast forcing
    plot(...
        time,H_C_fast,...
        'LineWidth',2,'Color',fainted_red)
    plot(...
        time(1:ind), H_C_fast(1:ind),...
        'LineWidth',2,'Color',red)
    plot(...
        time(ind),H_C_fast(ind),'.',...
        'MarkerSize',25,'Color',red)   
    box off
    xlim([-1000 3200])
    ylim([0 0.4])
    xticks([-1000 0 1000 2000 3000])
    xticklabels({})
    yticks([0 0.2 0.4])
    yticklabels([0.0 0.2 0.4])
    ylabel(['Freshwater ';'Hosing [Sv]']);

    % second subplot (salinity time series)
    h2 = subplot(7,13,[53:58,66:71,79:84]);
    hold on

    % QSE
    plot(...
        time,Q_on_fast,'--k',...
        time,Q_off_fast,'--k',...
        'LineWidth',1,'Color',[0.75 0.75 0.75])

    %time series 
    plot(...
        time,Q_slow,'Color',fainted_green,...
        'LineWidth',2)
    plot(...
        time(end),Q_slow(end),'.','Color', ...
        fainted_green,'MarkerSize',25)
    plot(...
        time(1:ind),Q_fast(1:ind),'Color',red,...
        'LineWidth',2)
    plot(...
        time(ind),Q_fast(ind),'.','Color',red,...
        'MarkerSize',25)

    xlim([-1000 3200])
    ylim([-10 15])
    xticks([-1000 0 1000 2000 3000])
    xticklabels([0 1000 2000 3000 4000])
    yticks([-10 0 10])
    yticklabels([-10 0 10])
    
    xlabel('time [years]')
    ylabel(['AMOC         ';'Strength [Sv]']);

    % third subplot (phase space)
    h3 = subplot(7,13,[8:13,21:26,34:39,60:65,73:78,86:91]);
    hold on
    
    % fast forcing
    if time(ind)<=0
        plot(...
            w_f(:,1),w_f(:,2),'-',...
            'LineWidth',2,'Color',0.75*[1 1 1])
    else
        plot(...
            w_f(500:end,1),w_f(500:end,2),'-',...
            'LineWidth',2,'Color',0.75*[1 1 1])
    end
    plot(...
        x_slow(:,1),x_slow(:,2),...
        'LineWidth',2,'Color',fainted_green)
    plot(...
        x_slow(end,1),x_slow(end,2),'.',...
        'Color',green,'MarkerSize',25)
    plot(...
        x_fast(1:ind,1),x_fast(1:ind,2),...
        'LineWidth',2,'Color',red)
    plot(...
        x_fast(1,1),x_fast(1,2),'.k',...
        'MarkerSize',30)
    plot(...
        x_fast(ind,1),x_fast(ind,2),'.',...
        'MarkerSize',25,'Color',red)

    plot(...
        S_on(1),S_on(2),'k.',...
        'MarkerSize',30)
    plot(...
        S_off(1),S_off(2),'k.',...
        'MarkerSize',30)
    plot(...
        S_sad(1),S_sad(2),'ko',... 
        'MarkerSize',6.5)

    xlabel('North Atlantic Salinity [psu]')
    ylabel('Tropical Atlantic Salinity [psu]')
    
    set(h1,'FontSize',20)
    set(h1,'Position',[0.09 0.59 0.43 0.33])

    set(h2,'FontSize',20)
    set(h2,'Position',[0.09 0.11 0.43 0.33])

    set(h3,'FontSize',20)
    set(h3,'Position',[0.62 0.11 0.35 0.82])


    xlim([-0.25 0.05]);
    ylim([0 0.3]);
    xticks([-0.2 -0.1 0])
    xticklabels([33 34 35])
    yticks([0 0.1 0.2 0.3])
    yticklabels([35 36 37 38])
    box on

    left   =  .3e3;
    bottom =  0.09e3;
    width  =  1.6e3;
    height =  .7e3;

    set(gcf,'position', [left bottom width height])
    set(gcf,'color','w');  
    set(gcf, 'renderer', 'painters')
    writeVideo(V_full,getframe(gcf));
end

% third part (fast, slow and critical)
for ind = 200:length(time)
    % the on- off- and saddle- states
    S_on  = [up_bran1(H_crit(time(ind)));...
        up_bran2(H_crit(time(ind)))];
    S_off = [lo_bran1(H_crit(time(ind)));...
        lo_bran2(H_crit(time(ind)))];
    S_sad =  [un_bran1(H_crit(time(ind)));...
        un_bran2(H_crit(time(ind)))];

    % tipping boundary -- the stable manifold on S_sad
    HH = @(tt)H_crit(time(ind));
    [~,w_c] = ode45 (@(t,x)BoxModel_2DH_IVP(t,x,HH), ...
        7000:-10:0,S_sad-[0; 0.001]);

    figure(10);
    clf
    annotation('textbox',...
        [0.412 0.94 0.266 0.06],...
        'String',{'Critical freshwater hosing'},...
        'LineStyle','none',...
        'Interpreter','latex',...
        'FontSize',20);

    % first subplot (forcing time series)
    h1 = subplot(7,13,[1:6,14:19,27:32]);
    hold on
    % slow forcing 
    plot(...
        time,H_C_slow,...
        'LineWidth',2,'Color',fainted_green)
    % fast forcing
    plot(...
        time,H_C_fast,...
        'LineWidth',2,'Color',fainted_red)  
    % critical forcing 
    plot(...
        time,H_C_crit,...
        'LineWidth',2,'Color',fainted_blue)
    plot(...
        time(1:ind), H_C_crit(1:ind),...
        'LineWidth',2,'Color',blue)
    plot(...
        time(ind),H_C_crit(ind),'.',...
        'MarkerSize',25,'Color',blue)

    box off
    xlim([-1000 3200])
    ylim([0 0.4])
    xticks([-1000 0 1000 2000 3000])
    xticklabels({})
    yticks([0 0.2 0.4])
    yticklabels([0.0 0.2 0.4])
    ylabel(['Freshwater ';'Hosing [Sv]']);

    % second subplot (salinity time series)
    h2 = subplot(7,13,[53:58,66:71,79:84]);
    hold on

    % QSE
    plot(...
        time,Q_on_crit,'--k',...
        time,Q_off_crit,'--k',...
        'LineWidth',1,'Color',[0.7 0.7 0.7])

    %time series 
    plot(...
        time,Q_slow,'Color',fainted_green,...
        'LineWidth',2)
    plot(...
        time(end),Q_slow(end),'.','Color', ...
        fainted_green,'MarkerSize',25)
    plot(...
        time,Q_fast,'Color',fainted_red,...
        'LineWidth',2)
    plot(...
        time(end),Q_fast(end),'.','Color', ...
        fainted_red,'MarkerSize',25)
    plot(...
        time(1:ind),Q_crit(1:ind),'Color',blue,...
        'LineWidth',2)
    plot(...
        time(ind),Q_crit(ind),'.','Color',blue,...
        'MarkerSize',25)

    xlim([-1000 3200])
    ylim([-10 15])
    xticks([-1000 0 1000 2000 3000])
    xticklabels([0 1000 2000 3000 4000])
    yticks([-10 0 10])
    yticklabels([-10 0 10])
    
    xlabel('time [years]')
    ylabel(['AMOC         ';'Strength [Sv]']);

    % third subplot (phase space)
    h3 = subplot(7,13,[8:13,21:26,34:39,60:65,73:78,86:91]);
    hold on
    
    % critical forcing

    if time(ind)<=-5
        plot(...
            w_c(:,1),w_c(:,2),'-','Color',0.75.*[1 1 1],...
            'LineWidth',2)
    else
        plot(...
            w_c(500:end,1),w_c(500:end,2),'-','Color',0.75.*[1 1 1],...
            'LineWidth',2)
    end

    plot(...
        x_slow(:,1),x_slow(:,2),'Color',fainted_green,...
        'LineWidth',2)
    plot(...
        x_slow(end,1),x_slow(end,2),'.','Color',fainted_green,...
        'MarkerSize',25)
    plot(...
        x_fast(end,1),x_fast(end,2),'.','Color',fainted_red,...
        'MarkerSize',25)
    plot(x_fast(:,1),x_fast(:,2),'Color',fainted_red,...
        'LineWidth',2)    
    plot(...
        x_crit(1:ind,1),x_crit(1:ind,2),'-',...
        'Color',blue,'LineWidth',2)
    plot(x_crit(1,1),x_crit(1,2),'.k',...
        'MarkerSize',30)
    plot(...
        x_crit(ind,1),x_crit(ind,2),'.',...
        'Color',blue,'MarkerSize',25)

    plot(...
        S_on(1),S_on(2),'k.',...
        'MarkerSize',30)
    plot(...
        S_off(1),S_off(2),'k.',...
        'MarkerSize',30)
    plot(...
        S_sad(1),S_sad(2),'ko',... 
        'MarkerSize',6.5)

    xlabel('North Atlantic Salinity [psu]')
    ylabel('Tropical Atlantic Salinity [psu]')
    
    set(h1,'FontSize',20)
    set(h1,'Position',[0.09 0.59 0.43 0.33])

    set(h2,'FontSize',20)
    set(h2,'Position',[0.09 0.11 0.43 0.33])

    set(h3,'FontSize',20)
    set(h3,'Position',[0.62 0.11 0.35 0.82])


    xlim([-0.25 0.05]);
    ylim([0 0.3]);
    xticks([-0.2 -0.1 0])
    xticklabels([33 34 35])
    yticks([0 0.1 0.2 0.3])
    yticklabels([35 36 37 38])
    box on

    left   =  .3e3;
    bottom =  0.09e3;
    width  =  1.6e3;
    height =  .7e3;

    set(gcf,'position', [left bottom width height])
    set(gcf,'color','w');  
    set(gcf, 'renderer', 'painters')
    writeVideo(V_full,getframe(gcf));
end
close(V_full)
%% test the critical rate 
% figure(1)
%   plot(...
%         x_crit(:,1),x_crit(:,2),'-',...
%         'Color','b','LineWidth',2)



