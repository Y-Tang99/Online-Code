% This code estimates the effective distance for ground motion simulations
% The basic idea is to make the left-hand side equals to the right-hand side 
% of the equation. Detailed information can be found in 
% "Boore, D. M. (2009). Comparing stochastic point-source and finite-source 
% ground-motion simulations: SMSIM and EXSIM. Bull.Seismol.Soc.Am.99:3202-3216."


Reff=53.2;      % this is final value you need to find


fQ=10;              % Boore (2009), change with G & Q


EL=REFF(Reff,fQ);    % Equation (6) Left-hand side in Boore (2009)

M=6.1;               % Moment magnitude

stress=50;        % Stress drop, bars
stress_ref=70;     % reference stress drop, used for finding FL & FW

%--------------------------------------------------------------------------
% Finite-fault Input Parameters (custom defined)
%--------------------------------------------------------------------------
 
FaultLat=27.11;            % latitude of upper edge of fault
FaultLon=103.35;            % longitude of upper edge of fault
Fstrike=162;                    % fault strike,degree (°)
Fdip=86;                      % fault dip, degree (°)
Rake=45;                      % rake angle, degree (°)

% Fault Dimensions 

 FL=42;       % fault length
 FW=20;       % fault width

% The following is used for the situation where the fault length and width
% are unknown. Wells and Coppersmith's correlation (1994) is used here

% if Rake == 0 || Rake == 180                                  %% Strike Slip
%     FL=10^(-2.57+0.62*M)*(stress_ref/stress)^(1/3);
%     FW=10^(-0.76+0.27*M)*(stress_ref/stress)^(1/3);
% elseif (Rake > 0) && (Rake < 180) && Fdip ~=0 && Fdip ~= 90  %% Reverse
%     FL=10^(-2.42+0.58*M)*(stress_ref/stress)^(1/3);
%     FW=10^(-1.61+0.41*M)*(stress_ref/stress)^(1/3);
% elseif (Rake > -180) && (Rake < 0) && Fdip ~=0 && Fdip ~= 90 %% Normal
%     FL=10^(-1.88+0.50*M)*(stress_ref/stress)^(1/3);
%     FW=10^(-1.14+0.35*M)*(stress_ref/stress)^(1/3);
% else                                                         %% Undefined
%     FL=10^(-2.44+0.59*M)*(stress_ref/stress)^(1/3);
%     FW=10^(-1.01+0.32*M)*(stress_ref/stress)^(1/3);
% end

% The following parameters are needed to locate the origin point

s1f=-FL/2;         % Along strike near edge            
s2f=FL/2;        % Along strike far edge 
w1f=-FW/2;         % Down dip near edge
w2f=FW/2;        % Down dip far edge
h_ref=10;    % fault depth to upper edge

h_min=3.0;     % Campbell depth to seismogenic region, usually set as 3.0

% Subfault Dimension

dl=2;                      % subfault length, no less than 1.5 km
dw=2;                      % subfault width, no less than 1.5 km

nl=round(FL/dl);           % number of subfaults along strike
nw=round(FW/dw);           % number of subfaults along dip
if nl<=1
    nl=1;
end
if nw<=1
    nw=1;
end

NF=nl*nw;                  % total number of subfaults

%--------------------------------------------------------------------------
% (3) Site Input Parameters
%--------------------------------------------------------------------------

% Site location

% Two options for determining site location: lattitude & longitude (LatLon),
% and distance & azimuth (DistAz). Users need to choose one for their purposes.

 SLIndex='LatLon';         % Latitude and longitude
% SLIndex='DistAz';          % Distance and Azimuth

%    SL1=5;                    % Input values to get site location
%    SL2=0;


% 51LZT
% SL1=28.899;              % site/station latitude
% SL2=105.4;        % site/station longitude

% 51YBY
% SL1=29;            % site/station latitude
% SL2=104.599;           % site/station longitude

% 53DTB
% SL1=26.399;        % site/station latitude
% SL2=103;             % site/station longitude

% 53DTD
% SL1=26.2;        % site/station latitude
% SL2=103.099;           % site/station longitude

% 53LDC
% SL1=27.2;           % site/station latitude
% SL2=103.599;        % site/station longitude

% 53QJT
% SL1=26.899;           % site/station latitude
% SL2=102.9;          % site/station longitude

% 53SFX
% SL1=28.6;           % site/station latitude
% SL2=104.4;          % site/station longitude

% 53ZTT
SL1=27.299;           % site/station latitude
SL2=103.699;          % site/station longitude



[SiteLat,SiteLon,R,Az]=FUNSL(SLIndex,SL1,SL2,FaultLat,FaultLon);

subR(nl,nw)=zeros();
ER0(nl,nw)=zeros();    % Equation (6) Right part in Boore (2009)

for i=1:1:nl
    for j=1:1:nw
        subR(i,j)=FUNsubR(R,h_ref,Fdip,Fstrike,Az,dl,dw,i,j);
        ER0(i,j)=REFF(subR(i,j),fQ);
    end
end

ER=sqrt(sum(sum(ER0.^2))/(NF));


function [Eq]=REFF(Reff,f)

% Geometric spreading function

% if Reff<=70
%     G=Reff^(-1.3);
% else
%     if Reff<=140
%         G=(70^(-0.2)/70^(1.3))*(Reff^0.2);
%     else
%         G=(70^(-0.2)/70^(1.3))*(140^(0.5)/140^(-0.2))*(Reff^(-0.5));
%     end
% end

if Reff<=50
    G=Reff^(-1.0);
elseif Reff<=90
    G=((50^0.3)/(50^1.0))*((Reff)^(-0.3));
elseif Reff<=120
    G=((50^0.3)/(50^1.0))*((90^1.1)/(90^0.3))*(Reff^(-1.1));
else
    G=((50^0.3)/(50^1.0))*((90^1.1)/(90^0.3))*((120^0.5)/(120^1.1))*(Reff^(-0.5));
end


cq=3.5;        % note this is different from beta0!!
Q0=180;
nq=0.5;
%Qmin=60;
%Q=max(Qmin,Q0*(subf.^(nq)));
Q=Q0*(f.^nq);

Ae1=-pi*f.*Reff;
Ae2=Q.*cq;
Ae=exp(Ae1./Ae2);
   
% cq=3.7;        % note this is different from beta0!!
% Q0=893;
% nq=0.32;
% Qmin=1000;
% Q=max(Qmin,Q0*(f.^(nq)));
% 
% Ae1=-pi*f.*Reff;
% Ae2=Q.*cq;
% Ae=exp(Ae1./Ae2);

Eq=G.*Ae;
end

function [SiteLat,SiteLon,R,Az] = FUNSL(SLIndex,SL1,SL2,FaultLat,FaultLon)
% This function is used to determine the site locations

if SLIndex == 'LatLon'
    SiteLat=SL1;
    SiteLon=SL2;
    % calculate the distance and azimuth of the site with rspect to the origin
    [ArcLen,~]=distance(SiteLat,SiteLon,FaultLat,FaultLon);
    R=ArcLen*6371*pi/180;    % epicentral distance
    R1=6371*(SiteLat-FaultLat)*pi/180;
    if R1>=R
        Az=180;
    else
        Az=acos(R1/R)*180/pi;
    end
    if SiteLon-FaultLon<=0
        Az=360-Az;
    end
else
    if SLIndex == 'DistAz'
        R=SL1;
        Az=SL2;
        ArcLen=(R/6371)*180/pi;
        [SiteLat,SiteLon]=reckon(FaultLat,FaultLon,ArcLen,Az);
    end
end

end

function [subR] = FUNsubR(R,h_Ref,Fdip,Fstrike,Az,dl,dw,i,j)
 % This function is used for finding the subfault distance
Fstrike_radians=(Az-Fstrike)*pi/180;
Fdip_radians=(90-Fdip)*pi/180;

t1=R*cos(Fstrike_radians)-(2*i-1)*dl/2;
t2=R*sin(Fstrike_radians)-((2*j-1)*dw/2)*sin(Fdip_radians);
t3=-h_Ref-((2*j-1)*dw/2)*cos(Fdip_radians);
subR=sqrt(t1^2+t2^2+t3^2);
end




