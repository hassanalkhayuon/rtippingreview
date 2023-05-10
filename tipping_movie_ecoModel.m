% A movie of the phase digrams of the plant-herbivore model (O'Keeffe and Wieczorek) 
% showing R-tipping with monotonic enviormantal conditions. i.e. increasing
% mortality rate for herbivore and growth rate for plants 
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
tspan_slow = -200:0.5:150;
tspan_fast = -200:0.5:150;
tspan_crit = -200:0.5:150;
 
rho_0 = 0.5;


r_slow = 0.03; %(track)
r_crit = 0.04315812; %critical), 
r_fast = 0.1; %(tipping).
drho  = 0.1;

rho_slow      = @(tt)((rho_0 + drho * sech(r_slow*tt))*(tt<0) + (rho_0 + drho)*(tt>=0));
rho_crit      = @(tt)((rho_0 + drho * sech(r_crit*tt))*(tt<0) + (rho_0 + drho)*(tt>=0));
rho_fast      = @(tt)((rho_0 + drho * sech(r_fast*tt))*(tt<0) + (rho_0 + drho)*(tt>=0));

init     = [Cs_bran1(rho_0),Cs_bran2(rho_0)]; 

odefun_slow       = @(tt,x)ecoODE(tt,x,rho_slow);
[~,x_slow]        = ode45(odefun_slow,tspan_slow,init,opts);

odefun_crit       = @(tt,x)ecoODE(tt,x,rho_crit);
[~,x_crit]        = ode45(odefun_crit,tspan_crit,init,opts);

odefun_fast       = @(tt,x)ecoODE(tt,x,rho_fast);
[time,x_fast]        = ode45(odefun_fast,tspan_fast,init,opts);
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
    rho_C_slow (ind)  = rho_slow(time(ind));
    H_QS_C_slow (ind)  = Cs_bran2(rho_slow(time(ind)));
    H_QS_E_slow (ind) =  H_bran2(rho_slow(time(ind)));

    rho_C_crit (ind)  = rho_crit(time(ind));
    H_QS_C_crit (ind)  = Cs_bran2(rho_crit(time(ind)));
    H_QS_E_crit (ind) =  H_bran2(rho_crit(time(ind)));


    rho_C_fast (ind)  = rho_fast(time(ind));
    H_QS_C_fast (ind)  = Cs_bran2(rho_fast(time(ind)));
    H_QS_E_fast (ind) =  H_bran2(rho_fast(time(ind)));
end


%% Slow forcing movie 

movie_name_s = 'PH_Rtipping_s';
V_slow = VideoWriter(movie_name_s,'MPEG-4');
V_slow.FrameRate = 30;
V_slow.Quality = 100;

open(V_slow);
end_of_run = length(time);
for ind = 1:end_of_run

    % the coexistence, plants-only and saddle state
    C     = [Cs_bran1(rho_slow(time(ind))); Cs_bran2(rho_slow(time(ind)))];
    P     = [H_bran1(rho_slow(time(ind))); H_bran2(rho_slow(time(ind)))];
    C_sad = [Cu_bran1(rho_slow(time(ind))); Cu_bran2(rho_slow(time(ind)))];

    % tipping boundary -- the stable manifold on S_sad
    Rho = @(tt)rho_slow(time(ind));
    [~,w_s1] = ode45 (@(t,x)ecoODE(t,x,Rho), ...
        500:-0.5:0,C_sad-[0.01; 0],opts);
    [~,w_s2] = ode45 (@(t,x)ecoODE(t,x,Rho), ...
        500:-0.5:0,C_sad+[0.01; 0],opts);
    

    figure(10);
    clf

    % first subplot (forcing time series)
    h1 = subplot(7,13,[1:6,14:19,27:32]);
    hold on

    % slow forcing 
    plot(...
        time,rho_C_slow,...
        'LineWidth',2,'Color',fainted_green)
    plot(...
        time(1:ind), rho_C_slow(1:ind),...
        'LineWidth',2,'Color',green)
    plot(...
        time(ind),rho_C_slow(ind),'.',...
        'MarkerSize',25,'Color',green)   
    box off
    xlim([-200 150])
    ylim([0.5 0.61])
    xticks([-200 -100 0 100 200])
    xticklabels({})
    yticks([0.5 0.55 0.6])
    yticklabels({'0.00' '0.25' '0.50'})
    ylabel(['Environmental';'conditions   ']);

    % second subplot (variable time series)
    h2 = subplot(7,13,[53:58,66:71,79:84]);
    hold on

    % QSE
    plot(...
        time,H_QS_C_slow,'--k',...
        time,H_QS_E_slow,'--k',...
        'LineWidth',1,'Color',[0.7 0.7 0.7])

    %time series 
    plot(...
        time(1:ind),x_slow(1:ind,2),'Color',green,...
        'LineWidth',2)
    plot(...
        time(ind),x_slow(ind,2),'.','Color',green,...
        'MarkerSize',25)

    xlim([-200 150])
    ylim([-1 11])
    xticks([-200 -100 0 100 200])
    xticklabels([0 100 200 300 400])
    yticks([0 3 6 9])
    yticklabels([0 3 6 9])
    
    xlabel('time [days]')
    ylabel(['Herbivore         ';'biomass [g/m$^2$] ']);

    % third subplot (phase space)
    h3 = subplot(7,13,[8:13,21:26,34:39,60:65,73:78,86:91]);
    hold on
    
    % slow forcing
    plot(...
        w_s1(:,1),w_s1(:,2),'-',...
        'LineWidth',2,'Color',0.75*[1 1 1])
    plot(...
        w_s2(:,1),w_s2(:,2),'-',...
        'LineWidth',2,'Color',0.75*[1 1 1])
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
        C(1),C(2),'k.',...
        'MarkerSize',30)
    plot(...
        P(1),P(2),'k.',...
        'MarkerSize',30)
    plot(...
        C_sad(1),C_sad(2),'ko',...
        'MarkerSize',6.5,'LineWidth',2)

    xlabel('Plant biomass [g/m$^2$]')
    ylabel('Herbivore biomass [g/m$^2$]')
    
    set(h1,'FontSize',20)
    set(h1,'Position',[0.09 0.6 0.43  0.33])

    set(h2,'FontSize',20)
    set(h2,'Position',[0.09 0.12 0.43 0.33])

    set(h3,'FontSize',20)
    set(h3,'Position',[0.62 0.12 0.35  0.82])


    xlim([8 32]);
    ylim([0 12]);
    xticks([10 15 20 25 30])
    xticklabels([10 15 20 25 30])
    yticks([0 3 6 9 12])
    yticklabels([0 3 6 9 12])
    box on

    left   =  .3e3;
    bottom =  0.09e3;
    width  =  1.6e3;
    height =  .7e3;

    set(gcf,'position', [left bottom width height])
    set(gcf,'color','w');  
    set(gcf, 'renderer', 'painters')
    writeVideo(V_slow,getframe(gcf));
end
close(V_slow);
%% Fast forcing movie 

movie_name_f = 'PH_Rtipping_f';
V_fast = VideoWriter(movie_name_f,'MPEG-4');
V_fast.FrameRate = 30;
V_fast.Quality = 100;

open(V_fast);
end_of_run = length(time);
for ind = 1:end_of_run
    % the coexistence, plants-only and saddle state
    C     = [Cs_bran1(rho_fast(time(ind))); Cs_bran2(rho_fast(time(ind)))];
    P     = [H_bran1(rho_fast(time(ind))); H_bran2(rho_fast(time(ind)))];
    C_sad = [Cu_bran1(rho_fast(time(ind))); Cu_bran2(rho_fast(time(ind)))];

    % tipping boundary -- the stable manifold on S_sad
    Rho = @(tt)rho_fast(time(ind));
    [~,w_f1] = ode45 (@(t,x)ecoODE(t,x,Rho), ...
        500:-0.5:0,C_sad-[0.01; 0],opts);
    [~,w_f2] = ode45 (@(t,x)ecoODE(t,x,Rho), ...
        500:-0.5:0,C_sad+[0.01; 0],opts);
    

    figure(10);
    clf

    % first subplot (forcing time series)
    h1 = subplot(7,13,[1:6,14:19,27:32]);
    hold on

    % fast forcing 
    plot(...
        time,rho_C_slow,...
        'LineWidth',2,'Color',fainted_green)
    plot(...
        time,rho_C_fast,...
        'LineWidth',2,'Color',fainted_red)
    plot(...
        time(1:ind), rho_C_fast(1:ind),...
        'LineWidth',2,'Color',red)
    plot(...
        time(ind),rho_C_fast(ind),'.',...
        'MarkerSize',25,'Color',red)   
    box off
    xlim([-200 150])
    ylim([0.5 0.61])
    xticks([-200 -100 0 100 200])
    xticklabels({})
    yticks([0.5 0.55 0.6])
    yticklabels({'0.00' '0.25' '0.50'})
    ylabel(['Environmental';'conditions   ']);

    % second subplot (salinity time series)
    h2 = subplot(7,13,[53:58,66:71,79:84]);
    hold on

    % QSE
    plot(...
        time,H_QS_C_fast,'--k',...
        time,H_QS_E_fast,'--k',...
        'LineWidth',1,'Color',[0.7 0.7 0.7])

    %time series 
    plot(...
        time,x_slow(:,2),'Color',fainted_green,...
        'LineWidth',2)
    plot(...
        time(end),x_slow(end,2),'.','Color',fainted_green,...
        'MarkerSize',25)
    plot(...
        time(1:ind),x_fast(1:ind,2),'Color',red,...
        'LineWidth',2)
    plot(...
        time(ind),x_fast(ind,2),'.','Color',red,...
        'MarkerSize',25)

    xlim([-200 150])
    ylim([-1 11])
    xticks([-200 -100 0 100 200])
    xticklabels([0 100 200 300 400])
    yticks([0 3 6 9])
    yticklabels([0 3 6 9])
    
    xlabel('time [days]')
    ylabel(['Herbivore         ';'biomass [g/m$^2$] ']);

    % third subplot (phase space)
    h3 = subplot(7,13,[8:13,21:26,34:39,60:65,73:78,86:91]);
    hold on
    
    % slow forcing
    plot(...
        w_f1(:,1),w_f1(:,2),'-',...
        'LineWidth',2,'Color',0.75*[1 1 1])
    plot(...
        w_f2(:,1),w_f2(:,2),'-',...
        'LineWidth',2,'Color',0.75*[1 1 1])
    plot(...
        x_slow(:,1),x_slow(:,2),...
        'LineWidth',2,'Color',fainted_green)
    plot(...
        x_slow(end,1),x_slow(end,2),'.',...
        'MarkerSize',25,'Color',fainted_green)
    plot(...
        x_fast(1:ind,1),x_fast(1:ind,2),...
        'LineWidth',2,'Color',red)
    plot(...
        x_fast(1,1),x_fast(1,2),'k.',...
        'MarkerSize',30)
    plot(...
        x_fast(ind,1),x_fast(ind,2),'.',...
        'MarkerSize',25,'Color',red)

    plot(...
        C(1),C(2),'k.',...
        'MarkerSize',30)
    plot(...
        P(1),P(2),'k.',...
        'MarkerSize',30)
    plot(...
        C_sad(1),C_sad(2),'ko',...
        'MarkerSize',6.5,'LineWidth',2)

    xlabel('Plant biomass [g/m$^2$]')
    ylabel('Herbivore biomass [g/m$^2$]')
    
    set(h1,'FontSize',20)
    set(h1,'Position',[0.09 0.6 0.43  0.33])

    set(h2,'FontSize',20)
    set(h2,'Position',[0.09 0.12 0.43 0.33])

    set(h3,'FontSize',20)
    set(h3,'Position',[0.62 0.12 0.35  0.82])


    xlim([8 32]);
    ylim([0 12]);
    xticks([10 15 20 25 30])
    xticklabels([10 15 20 25 30])
    yticks([0 3 6 9 12])
    yticklabels([0 3 6 9 12])
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
close(V_fast);

%% Critical forcing movie 

movie_name_c = 'PH_Rtipping_c';
V_crit = VideoWriter(movie_name_c,'MPEG-4');
V_crit.FrameRate = 30;
V_crit.Quality = 100;

open(V_crit);
end_of_run = length(time);
for ind = 1:end_of_run
    % the coexistence, plants-only and saddle state
    C     = [Cs_bran1(rho_crit(time(ind))); Cs_bran2(rho_crit(time(ind)))];
    P     = [H_bran1(rho_crit(time(ind))); H_bran2(rho_crit(time(ind)))];
    C_sad = [Cu_bran1(rho_crit(time(ind))); Cu_bran2(rho_crit(time(ind)))];

    % tipping boundary -- the stable manifold on S_sad
    Rho = @(tt)rho_crit(time(ind));
    [~,w_c1] = ode45 (@(t,x)ecoODE(t,x,Rho), ...
        500:-0.5:0,C_sad-[0.01; 0],opts);
    [~,w_c2] = ode45 (@(t,x)ecoODE(t,x,Rho), ...
        500:-0.5:0,C_sad+[0.01; 0],opts);
    

    figure(10);
    clf

    % first subplot (forcing time series)
    h1 = subplot(7,13,[1:6,14:19,27:32]);
    hold on

    % critical forcing 
    plot(...
        time,rho_C_slow,...
        'LineWidth',2,'Color',fainted_green)
    plot(...
        time,rho_C_fast,...
        'LineWidth',2,'Color',fainted_red)
    plot(...
        time,rho_C_crit,...
        'LineWidth',2,'Color',fainted_blue)
    plot(...
        time(1:ind), rho_C_crit(1:ind),...
        'LineWidth',2,'Color',blue)
    plot(...
        time(ind),rho_C_crit(ind),'.',...
        'MarkerSize',25,'Color',blue)   
    box off
    xlim([-200 150])
    ylim([0.5 0.61])
    xticks([-200 -100 0 100 200])
    xticklabels({})
    yticks([0.5 0.55 0.6])
    yticklabels({'0.00' '0.25' '0.50'})
    ylabel(['Environmental';'conditions   ']);

    % second subplot (salinity time series)
    h2 = subplot(7,13,[53:58,66:71,79:84]);
    hold on

    % QSE
    plot(...
        time,H_QS_C_crit,'--k',...
        time,H_QS_E_crit,'--k',...
        'LineWidth',1,'Color',[0.7 0.7 0.7])

    %time series 
    plot(...
        time,x_slow(:,2),'Color',fainted_green,...
        'LineWidth',2)
    plot(...
        time(end),x_slow(end,2),'.','Color',fainted_green,...
        'MarkerSize',25)
    plot(...
        time,x_fast(:,2),'Color',fainted_red,...
        'LineWidth',2)
    plot(...
        time(end),x_fast(end,2),'.','Color',fainted_red,...
        'MarkerSize',25)
    plot(...
        time(1:ind),x_crit(1:ind,2),'Color',blue,...
        'LineWidth',2)
    plot(...
        time(ind),x_crit(ind,2),'.','Color',blue,...
        'MarkerSize',25)

    xlim([-200 150])
    ylim([-1 11])
    xticks([-200 -100 0 100 200])
    xticklabels([0 100 200 300 400])
    yticks([0 3 6 9])
    yticklabels([0 3 6 9])
    
    xlabel('time [days]')
    ylabel(['Herbivore         ';'biomass [g/m$^2$] ']);

    % third subplot (phase space)
    h3 = subplot(7,13,[8:13,21:26,34:39,60:65,73:78,86:91]);
    hold on
    
    % slow forcing
    plot(...
        w_c1(:,1),w_c1(:,2),'-',...
        'LineWidth',2,'Color',0.75*[1 1 1])
    plot(...
        w_c2(:,1),w_c2(:,2),'-',...
        'LineWidth',2,'Color',0.75*[1 1 1])
    plot(...
        x_slow(:,1),x_slow(:,2),...
        'LineWidth',2,'Color',fainted_green)
    plot(...
        x_slow(end,1),x_slow(end,2),'.',...
        'MarkerSize',25,'Color',fainted_green)
    plot(...
        x_fast(:,1),x_fast(:,2),...
        'LineWidth',2,'Color',fainted_red)
    plot(...
        x_fast(end,1),x_fast(end,2),'.',...
        'MarkerSize',25,'Color',fainted_red)
    plot(...
        x_crit(1:ind,1),x_crit(1:ind,2),...
        'LineWidth',2,'Color',blue)
    plot(...
        x_crit(1,1),x_crit(1,2),'k.',...
        'MarkerSize',30)
    plot(...
        x_crit(ind,1),x_crit(ind,2),'.',...
        'MarkerSize',25,'Color',blue)

    plot(...
        C(1),C(2),'k.',...
        'MarkerSize',30)
    plot(...
        P(1),P(2),'k.',...
        'MarkerSize',30)
    plot(...
        C_sad(1),C_sad(2),'ko',...
        'MarkerSize',6.5,'LineWidth',2)

    xlabel('Plant biomass [g/m$^2$]')
    ylabel('Herbivore biomass [g/m$^2$]')
    
    set(h1,'FontSize',20)
    set(h1,'Position',[0.09 0.6 0.43  0.33])

    set(h2,'FontSize',20)
    set(h2,'Position',[0.09 0.12 0.43 0.33])

    set(h3,'FontSize',20)
    set(h3,'Position',[0.62 0.12 0.35  0.82])


    xlim([8 32]);
    ylim([0 12]);
    xticks([10 15 20 25 30])
    xticklabels([10 15 20 25 30])
    yticks([0 3 6 9 12])
    yticklabels([0 3 6 9 12])
    box on

    left   =  .3e3;
    bottom =  0.09e3;
    width  =  1.6e3;
    height =  .7e3;

    set(gcf,'position', [left bottom width height])
    set(gcf,'color','w');  
    set(gcf, 'renderer', 'painters')
    writeVideo(V_crit,getframe(gcf));
end
close(V_crit);
%% the full movie 
movie_name_full = 'PH_Rtipping_full';
V_full = VideoWriter(movie_name_full,'MPEG-4');
V_full.FrameRate = 30;
V_full.Quality = 100;

open(V_full);

% slow forcing part
end_of_run = length(time);
for ind = 1:end_of_run

    % the coexistence, plants-only and saddle state
    C     = [Cs_bran1(rho_slow(time(ind))); Cs_bran2(rho_slow(time(ind)))];
    P     = [H_bran1(rho_slow(time(ind))); H_bran2(rho_slow(time(ind)))];
    C_sad = [Cu_bran1(rho_slow(time(ind))); Cu_bran2(rho_slow(time(ind)))];

    % tipping boundary -- the stable manifold on S_sad
    Rho = @(tt)rho_slow(time(ind));
    [~,w_s1] = ode45 (@(t,x)ecoODE(t,x,Rho), ...
        500:-0.5:0,C_sad-[0.01; 0],opts);
    [~,w_s2] = ode45 (@(t,x)ecoODE(t,x,Rho), ...
        500:-0.5:0,C_sad+[0.01; 0],opts);
    

    figure(10);
    clf
    annotation('textbox',...
        [0.348125 0.94 0.309875 0.06],...
        'String',{'Slow change in environmental conditions'},...
        'LineStyle','none',...
        'Interpreter','latex',...
        'FontSize',20);

    % first subplot (forcing time series)
    h1 = subplot(7,13,[1:6,14:19,27:32]);
    hold on

    % slow forcing 
    plot(...
        time,rho_C_slow,...
        'LineWidth',2,'Color',fainted_green)
    plot(...
        time(1:ind), rho_C_slow(1:ind),...
        'LineWidth',2,'Color',green)
    plot(...
        time(ind),rho_C_slow(ind),'.',...
        'MarkerSize',25,'Color',green)   
    box off
    xlim([-200 150])
    ylim([0.5 0.61])
    xticks([-200 -100 0 100 200])
    xticklabels({})
    yticks([0.5 0.55 0.6])
    yticklabels({'0.00' '0.25' '0.50'})
    ylabel(['Environmental';'conditions   ']);

    % second subplot (salinity time series)
    h2 = subplot(7,13,[53:58,66:71,79:84]);
    hold on

    % QSE
    plot(...
        time,H_QS_C_slow,'--k',...
        time,H_QS_E_slow,'--k',...
        'LineWidth',1,'Color',[0.7 0.7 0.7])

    %time series 
    plot(...
        time(1:ind),x_slow(1:ind,2),'Color',green,...
        'LineWidth',2)
    plot(...
        time(ind),x_slow(ind,2),'.','Color',green,...
        'MarkerSize',25)

    xlim([-200 150])
    ylim([-1 11])
    xticks([-200 -100 0 100 200])
    xticklabels([0 100 200 300 400])
    yticks([0 3 6 9])
    yticklabels([0 3 6 9])
    
    xlabel('time [days]')
    ylabel(['Herbivore         ';'biomass [g/m$^2$] ']);

    % third subplot (phase space)
    h3 = subplot(7,13,[8:13,21:26,34:39,60:65,73:78,86:91]);
    hold on
    
    % slow forcing
    plot(...
        w_s1(:,1),w_s1(:,2),'-',...
        'LineWidth',2,'Color',0.75*[1 1 1])
    plot(...
        w_s2(:,1),w_s2(:,2),'-',...
        'LineWidth',2,'Color',0.75*[1 1 1])
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
        C(1),C(2),'k.',...
        'MarkerSize',30)
    plot(...
        P(1),P(2),'k.',...
        'MarkerSize',30)
    plot(...
        C_sad(1),C_sad(2),'ko',...
        'MarkerSize',6.5,'LineWidth',2)

    xlabel('Plant biomass [g/m$^2$]')
    ylabel('Herbivore biomass [g/m$^2$]')
    
    set(h1,'FontSize',20)
    set(h1,'Position',[0.09 0.59 0.43 0.33])

    set(h2,'FontSize',20)
    set(h2,'Position',[0.09 0.11 0.43 0.33])

    set(h3,'FontSize',20)
    set(h3,'Position',[0.62 0.11 0.35 0.82])


    xlim([8 32]);
    ylim([0 12]);
    xticks([10 15 20 25 30])
    xticklabels([10 15 20 25 30])
    yticks([0 3 6 9 12])
    yticklabels([0 3 6 9 12])
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

% fast forcing part 

for ind = 1:end_of_run
    % the coexistence, plants-only and saddle state
    C     = [Cs_bran1(rho_fast(time(ind))); Cs_bran2(rho_fast(time(ind)))];
    P     = [H_bran1(rho_fast(time(ind))); H_bran2(rho_fast(time(ind)))];
    C_sad = [Cu_bran1(rho_fast(time(ind))); Cu_bran2(rho_fast(time(ind)))];

    % tipping boundary -- the stable manifold on S_sad
    Rho = @(tt)rho_fast(time(ind));
    [~,w_f1] = ode45 (@(t,x)ecoODE(t,x,Rho), ...
        500:-0.5:0,C_sad-[0.01; 0],opts);
    [~,w_f2] = ode45 (@(t,x)ecoODE(t,x,Rho), ...
        500:-0.5:0,C_sad+[0.01; 0],opts);
    

    figure(10);
    clf
    annotation('textbox',...
        [0.348125 0.94 0.309875 0.06],...
        'String',{'Fast change in environmental conditions'},...
        'LineStyle','none',...
        'Interpreter','latex',...
        'FontSize',20);

    % first subplot (forcing time series)
    h1 = subplot(7,13,[1:6,14:19,27:32]);
    hold on

    % fast forcing 
    plot(...
        time,rho_C_slow,...
        'LineWidth',2,'Color',fainted_green)
    plot(...
        time,rho_C_fast,...
        'LineWidth',2,'Color',fainted_red)
    plot(...
        time(1:ind), rho_C_fast(1:ind),...
        'LineWidth',2,'Color',red)
    plot(...
        time(ind),rho_C_fast(ind),'.',...
        'MarkerSize',25,'Color',red)   
    box off
    xlim([-200 150])
    ylim([0.5 0.61])
    xticks([-200 -100 0 100 200])
    xticklabels({})
    yticks([0.5 0.55 0.6])
    yticklabels({'0.0' '0.5' '1.0'})
    ylabel(['Environmental';'conditions   ']);

    % second subplot (salinity time series)
    h2 = subplot(7,13,[53:58,66:71,79:84]);
    hold on

    % QSE
    plot(...
        time,H_QS_C_fast,'--k',...
        time,H_QS_E_fast,'--k',...
        'LineWidth',1,'Color',[0.7 0.7 0.7])

    %time series 
    plot(...
        time,x_slow(:,2),'Color',fainted_green,...
        'LineWidth',2)
    plot(...
        time(end),x_slow(end,2),'.','Color',fainted_green,...
        'MarkerSize',25)
    plot(...
        time(1:ind),x_fast(1:ind,2),'Color',red,...
        'LineWidth',2)
    plot(...
        time(ind),x_fast(ind,2),'.','Color',red,...
        'MarkerSize',25)

    xlim([-200 150])
    ylim([-1 11])
    xticks([-200 -100 0 100 200])
    xticklabels([0 100 200 300 400])
    yticks([0 3 6 9])
    yticklabels([0 3 6 9])
    
    xlabel('time [days]')
    ylabel(['Herbivore         ';'biomass [g/m$^2$] ']);

    % third subplot (phase space)
    h3 = subplot(7,13,[8:13,21:26,34:39,60:65,73:78,86:91]);
    hold on
    
    % slow forcing
    plot(...
        w_f1(:,1),w_f1(:,2),'-',...
        'LineWidth',2,'Color',0.75*[1 1 1])
    plot(...
        w_f2(:,1),w_f2(:,2),'-',...
        'LineWidth',2,'Color',0.75*[1 1 1])
    plot(...
        x_slow(:,1),x_slow(:,2),...
        'LineWidth',2,'Color',fainted_green)
    plot(...
        x_slow(end,1),x_slow(end,2),'.',...
        'MarkerSize',25,'Color',fainted_green)
    plot(...
        x_fast(1:ind,1),x_fast(1:ind,2),...
        'LineWidth',2,'Color',red)
    plot(...
        x_fast(1,1),x_fast(1,2),'k.',...
        'MarkerSize',30)
    plot(...
        x_fast(ind,1),x_fast(ind,2),'.',...
        'MarkerSize',25,'Color',red)

    plot(...
        C(1),C(2),'k.',...
        'MarkerSize',30)
    plot(...
        P(1),P(2),'k.',...
        'MarkerSize',30)
    plot(...
        C_sad(1),C_sad(2),'ko',...
        'MarkerSize',6.5,'LineWidth',2)

    xlabel('Plant biomass [g/m$^2$]')
    ylabel('Herbivore biomass [g/m$^2$]')
    
    set(h1,'FontSize',20)
    set(h1,'Position',[0.09 0.59 0.43 0.33])

    set(h2,'FontSize',20)
    set(h2,'Position',[0.09 0.11 0.43 0.33])

    set(h3,'FontSize',20)
    set(h3,'Position',[0.62 0.11 0.35 0.82])


    xlim([8 32]);
    ylim([0 12]);
    xticks([10 15 20 25 30])
    xticklabels([10 15 20 25 30])
    yticks([0 3 6 9 12])
    yticklabels([0 3 6 9 12])
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

% critical forcing part

for ind = 1:end_of_run
    % the coexistence, plants-only and saddle state
    C     = [Cs_bran1(rho_crit(time(ind))); Cs_bran2(rho_crit(time(ind)))];
    P     = [H_bran1(rho_crit(time(ind))); H_bran2(rho_crit(time(ind)))];
    C_sad = [Cu_bran1(rho_crit(time(ind))); Cu_bran2(rho_crit(time(ind)))];

    % tipping boundary -- the stable manifold on S_sad
    Rho = @(tt)rho_crit(time(ind));
    [~,w_c1] = ode45 (@(t,x)ecoODE(t,x,Rho), ...
        500:-0.5:0,C_sad-[0.01; 0],opts);
    [~,w_c2] = ode45 (@(t,x)ecoODE(t,x,Rho), ...
        500:-0.5:0,C_sad+[0.01; 0],opts);
    

    figure(10);
    clf
    annotation('textbox',...
        [0.339375 0.94 0.328125 0.06],...
        'String',{'Critical change in environmental conditions'},...
        'LineStyle','none',...
        'Interpreter','latex',...
        'FontSize',20);


    % first subplot (forcing time series)
    h1 = subplot(7,13,[1:6,14:19,27:32]);
    hold on

    % critical forcing 
    plot(...
        time,rho_C_slow,...
        'LineWidth',2,'Color',fainted_green)
    plot(...
        time,rho_C_fast,...
        'LineWidth',2,'Color',fainted_red)
    plot(...
        time,rho_C_crit,...
        'LineWidth',2,'Color',fainted_blue)
    plot(...
        time(1:ind), rho_C_crit(1:ind),...
        'LineWidth',2,'Color',blue)
    plot(...
        time(ind),rho_C_crit(ind),'.',...
        'MarkerSize',25,'Color',blue)   
    box off
    xlim([-200 150])
    ylim([0.5 0.61])
    xticks([-200 -100 0 100 200])
    xticklabels({})
    yticks([0.5 0.55 0.6])
    yticklabels({'0.0' '0.5' '1.0'})
    ylabel(['Environmental';'conditions   ']);

    % second subplot (salinity time series)
    h2 = subplot(7,13,[53:58,66:71,79:84]);
    hold on

    % QSE
    plot(...
        time,H_QS_C_crit,'--k',...
        time,H_QS_E_crit,'--k',...
        'LineWidth',1,'Color',[0.7 0.7 0.7])

    %time series 
    plot(...
        time,x_slow(:,2),'Color',fainted_green,...
        'LineWidth',2)
    plot(...
        time(end),x_slow(end,2),'.','Color',fainted_green,...
        'MarkerSize',25)
    plot(...
        time,x_fast(:,2),'Color',fainted_red,...
        'LineWidth',2)
    plot(...
        time(end),x_fast(end,2),'.','Color',fainted_red,...
        'MarkerSize',25)
    plot(...
        time(1:ind),x_crit(1:ind,2),'Color',blue,...
        'LineWidth',2)
    plot(...
        time(ind),x_crit(ind,2),'.','Color',blue,...
        'MarkerSize',25)

    xlim([-200 150])
    ylim([-1 11])
    xticks([-200 -100 0 100 200])
    xticklabels([0 100 200 300 400])
    yticks([0 3 6 9])
    yticklabels([0 3 6 9])
    
    xlabel('time [days]')
    ylabel(['Herbivore         ';'biomass [g/m$^2$] ']);

    % third subplot (phase space)
    h3 = subplot(7,13,[8:13,21:26,34:39,60:65,73:78,86:91]);
    hold on
    
    % slow forcing
    plot(...
        w_c1(:,1),w_c1(:,2),'-',...
        'LineWidth',2,'Color',0.75*[1 1 1])
    plot(...
        w_c2(:,1),w_c2(:,2),'-',...
        'LineWidth',2,'Color',0.75*[1 1 1])
    plot(...
        x_slow(:,1),x_slow(:,2),...
        'LineWidth',2,'Color',fainted_green)
    plot(...
        x_slow(end,1),x_slow(end,2),'.',...
        'MarkerSize',25,'Color',fainted_green)
    plot(...
        x_fast(:,1),x_fast(:,2),...
        'LineWidth',2,'Color',fainted_red)
    plot(...
        x_fast(end,1),x_fast(end,2),'.',...
        'MarkerSize',25,'Color',fainted_red)
    plot(...
        x_crit(1:ind,1),x_crit(1:ind,2),...
        'LineWidth',2,'Color',blue)
    plot(...
        x_crit(1,1),x_crit(1,2),'k.',...
        'MarkerSize',30)
    plot(...
        x_crit(ind,1),x_crit(ind,2),'.',...
        'MarkerSize',25,'Color',blue)

    plot(...
        C(1),C(2),'k.',...
        'MarkerSize',30)
    plot(...
        P(1),P(2),'k.',...
        'MarkerSize',30)
    plot(...
        C_sad(1),C_sad(2),'ko',...
        'MarkerSize',6.5,'LineWidth',2)

    xlabel('Plant biomass [g/m$^2$]')
    ylabel('Herbivore biomass [g/m$^2$]')
    
    set(h1,'FontSize',20)
    set(h1,'Position',[0.09 0.59 0.43 0.33])

    set(h2,'FontSize',20)
    set(h2,'Position',[0.09 0.11 0.43 0.33])

    set(h3,'FontSize',20)
    set(h3,'Position',[0.62 0.11 0.35 0.82])

    xlim([8 32]);
    ylim([0 12]);
    xticks([10 15 20 25 30])
    xticklabels([10 15 20 25 30])
    yticks([0 3 6 9 12])
    yticklabels([0 3 6 9 12])
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

close(V_full);

%% test the critical rate 
% figure(1)
% plot(...
%     x_crit(:,1),x_crit(:,2),'-',...
%     'Color','b','LineWidth',2)
% C_sad = [Cu_bran1(rho_0+drho), Cu_bran2(rho_0+drho)];
% norm(x_crit(end,:) - C_sad)

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