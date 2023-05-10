function [ z ] = BoxModel_2DH_IVP(t,x,HH)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

S0 = 0.035; % dimensionless

alpha = 0.12; %kg m^-3 C^-1
beta = 790.0;%kg m^-3 
Y = 100*3.15e7; %sec/year

% H0 = par(1);
% r  = par(2);
% dH = par(3);


H = HH(t);%H0 + dH * sech(r*t);

% 1 x CO2
% VN = 0.3261e17;%m^3
% VT = 0.7777e17; %m^3
% VS = 0.8897e17; %m^3
% VIP = 2.2020e17; %m^3
% VB = 8.6490e17; %m^3
% 
% % SN = 0.034912;
% % ST = 0.035435;
% 
% Ts = 4.773; %C
% T0 = 2.65; %C
% 
% C = 4.4463e16; %m^3 (total sum of V*S)
% 
% lambda = 2.79e7; %m^6kg^-1s^-1
% gamma = 0.39; %dimensionless
% mu = 5.5e-8; %deg s m^-3
% 
% KN = 5.456e6; %m^3s^-1
% KS = 5.447e6; %m^3s^-1
% 
% FN = (0.384+0.070.*H).*1e6;
% FT = (-0.723+0.752.*H).*1e6; %m^3s^-1
%====================================================

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
gamma = 0.36;
mu = 22e-8; %deg s m^-3

KN = 1.762e6; %m^3s^-1
KS = 1.872e6; %m^3s^-1

FN = (0.486+0.1311.*H).*1e6;
FT = (-0.997+0.6961.*H).*1e6; %m^3s^-1
%====================================================

SN = x(1,:);
ST = x(2,:);
SS = (0.034427-S0).*100;
%SIP = (0.034668-S0).*100;
SB = (0.034538-S0).*100;
SIP = 100.*(C-(VN.*SN+VT.*ST+VS.*SS+VB.*SB)./100-S0.*(VB+VN+VT+VIP+VS))./VIP;

q = lambda.*(alpha.*(Ts-T0)+beta.*(x(1,:)./100-SS./100))/(1+lambda*alpha*mu);
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