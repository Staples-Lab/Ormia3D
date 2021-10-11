function dx=TympanaCoupling(t,x,ft,f,dt,theta,ModeFlag)
 ff = interp1(ft,f,t); % Interpolate the data set (ft,f) at time t
%This is the current function that models the tympanal response.
%dt=travel time between tympana (in air)
%t=timestep
%x=displacement
%ft=incoming sound wave signals (time vector)
%f=incoming sound wave signals (pressure vector)
%ModeFlag=type of simulation (Miles or q2D)
%Angle of origin of sound relative to fly midline

%Sound wave interpolation
 tf=isnan(ff);
 if tf == 1
     ff1=0;
     ff2=0;
 else
   if dt<0
     ff2=ff;
     ff1=interp1(ft,f,t-abs(dt));
     tf2=isnan(ff1);
     if tf2 == 1
       ff1 = 0;
     end
   elseif dt>0
     ff1=ff;
     ff2=interp1(ft,f,t-abs(dt));
     tf2=isnan(ff2);
     if tf2 == 1
       ff2 = 0;
     end
   elseif dt==0
     ff1=ff;
     ff2=ff;
     tf2=isnan(ff2);
     if tf2 == 1
       ff2 = 0;
       ff1=0;
     end
   end
 end

 %% Constant definitions for system
 AngleBreak=30;                           %Point of new implementation
 SaturationPoint=55;                      %Complete saturation point
 baseDash=1.15e-5;                        %Base dashpot value (from animal)
 BaseSpring=0.576;                        %Base spring value (from animal)
 maxSpringFactor=BaseSpring*3.5747;       %Saturated spring factor, established experimentally
 maxDashFactor=baseDash*6.5;              %Satured dashpot factor, established experimentally
 m=2.88e-10;                              %kg, total mass of moving parts, Miles et al. symmetric value
 s=0.288e-6;                              %m^2, Tympanal area (1 tympanum)
 s2=s;
 s1=s;
 f1=ff1*s1;                               %Newtons, max. force on right membrane, SCS (larger)
 f2=ff2*s2;                               %Newtons, max. force on left membrane, SCS (smaller)
 %% Slope equations
 % Spring
 slopeValSpring=(maxSpringFactor-BaseSpring)/(SaturationPoint-AngleBreak);
 OffSetBSpring=BaseSpring-slopeValSpring*AngleBreak;
 Springval=slopeValSpring*abs(theta)+OffSetBSpring;
 %Damper
 slopeValDash=(maxDashFactor-baseDash)/(SaturationPoint-AngleBreak);
 OffSetBDash=baseDash-slopeValDash*AngleBreak;
 DashVal=slopeValDash*abs(theta)+OffSetBDash;


%% End updated section for model / novel behavior (q2D)
if ModeFlag==1
%% Begin section for determining model implementation (values)
 if abs(theta)<=AngleBreak %Interior to 'switch' (normal behavior)
     k1=BaseSpring; % N/m, right spring constant, Miles et al. symmetric value
     k2=k1; % N/m, left spring constant, Miles et al. symmetric value

     c1=baseDash; % N-s/m, right dashpot constant, Miles et al. symmetric value
     c2=c1; % N-s/m, left dashpot constant, Miles et al. symmetric value
     k3=5.18; % N/m, pivot spring constant, Miles et al. symmetric value
     c3=2.88e-5; % N-s/m, pivot dashpot constant, Miles et al. symmetric value
%     left side negative response ramp
elseif theta>AngleBreak && theta<SaturationPoint
     k2=Springval; % N/m, right spring constant, Miles et al. symmetric value
     k1=BaseSpring; % N/m, left spring constant, Miles et al. symmetric value
     c2=DashVal; % N-s/m, right dashpot constant, Miles et al. symmetric value
     c1=baseDash; % N-s/m, left dashpot constant, Miles et al. symmetric value
     k3=5.18; % N/m, pivot spring constant, Miles et al. symmetric value
     c3=2.88e-5; % N-s/m, pivot dashpot constant, Miles et al. symmetric value
%     right side negative response ramp
 elseif theta<-AngleBreak && theta>-SaturationPoint
     k2=BaseSpring; % N/m, right spring constant, Miles et al. symmetric value
     k1=Springval; % N/m, left spring constant, Miles et al. symmetric value
     c2=baseDash; % N-s/m, right dashpot constant, Miles et al. symmetric value
     c1=DashVal; % N-s/m, left dashpot constant, Miles et al. symmetric value
     k3=5.18; % N/m, pivot spring constant, Miles et al. symmetric value
     c3=2.88e-5; % N-s/m, pivot dashpot constant, Miles et al. symmetric value
 elseif theta>=SaturationPoint
     k2=maxSpringFactor; % N/m, right spring constant, Miles et al. symmetric value
     k1=BaseSpring; % N/m, left spring constant, Miles et al. symmetric value
     c2=maxDashFactor; % N-s/m, right dashpot constant, Miles et al. symmetric value
     c1=baseDash; % N-s/m, left dashpot constant, Miles et al. symmetric value
     k3=5.18; % N/m, pivot spring constant, Miles et al. symmetric value
     c3=2.88e-5; % N-s/m, pivot dashpot constant, Miles et al. symmetric value
 elseif theta<=-SaturationPoint
     k2=BaseSpring; % N/m, right spring constant, Miles et al. symmetric value
     k1=maxSpringFactor; % N/m, left spring constant, Miles et al. symmetric value
     c2=baseDash; % N-s/m, right dashpot constant, Miles et al. symmetric value
     c1=maxDashFactor; % N-s/m, left dashpot constant, Miles et al. symmetric value
     k3=5.18; % N/m, pivot spring constant, Miles et al. symmetric value
     c3=2.88e-5; % N-s/m, pivot dashpot constant, Miles et al. symmetric value
 end

elseif ModeFlag==0
 k1=BaseSpring; % N/m, right spring constant, Miles et al. symmetric value
 k2=BaseSpring; % N/m, left spring constant, Miles et al. symmetric value
 c1=baseDash; % N-s/m, right dashpot constant, Miles et al. symmetric value
 c2=baseDash; % N-s/m, left dashpot constant, Miles et al. symmetric value
 k3=5.18; % N/m, pivot spring constant, Miles et al. symmetric value
 c3=2.88e-5; % N-s/m, pivot dashpot constant, Miles et al. symmetric value
end


 %% ODE Model
 dx(1)=x(3);
 dx(2)=x(4);
 dx(3)=(1/m)*(f1-(k1+k3)*x(1)-k3*x(2)-(c1+c3)*x(3)-c3*x(4));  %RIGHT MEMBRANE (fly's POV) 1= right
 dx(4)=(1/m)*(f2-k3*x(1)-(k2+k3)*x(2)-c3*x(3)-(c2+c3)*x(4));  %LEFT MEMBRANE (fly's POV) 2= left
 dx=dx';
end
