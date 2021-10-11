clc; clear variables; close all;
ModeFlag=1;                                 %ModeFlag determines whether or not we're using the adaptive mode. 1= new q2D mode, 0 is standard Miles mode
FreqVal=6000;                               %Incident sound frequency
d=0.0012;                                   %Distance between points
AngleVal=45;                                %Incident sound angle
dt=(d/344)*sind(AngleVal);                  %Natural signal delay
ft=linspace(0,0.005,2000);                  %Time vector for sound
f=72*sin(2*pi*FreqVal*ft);                  %Incident sound

%% ODE Simulation Section
    x0=[0 0 0 0];                           %Initial conditions for ode simulation
    tspan=[0 0.005];                        %timespan for ode simulation
    [t,x]=ode45(@(t,x)TympanaCoupling(t,x,ft,f,dt,AngleVal,ModeFlag),tspan,x0);
