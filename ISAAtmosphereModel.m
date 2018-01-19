function [T,a,P,rho,mu] = ISAAtmosphereModel( h, V, ~ )
%ISAAtmosphereModel determines the standard T,rho,P,a,mu or mach number
%based on the standard altitude, h [m] OR [ft]
%
%   [T,a,P,rho]=ISAAatmosphereModel(h) 
%   returns T [Kelvin], a [m/s], P [Pa],and rho [kg/m^3]

%   [T,a,P,rho]=ISAAatmosphereModel(h,1) 
%   returns T [Rankine, a [ft/s], P [lbf/ft^2], and rho [lbf/ft^3]

%   [M]=ISAAatmosphereModel(h,v)
%   returns Mach number v/a at altitude, h, where a is in units m/s

%   [M]=ISAAatmosphereModel(h,v,1)
%   returns Mach number v/a at altitude, h, where a is in units ft/s

%   [T,a,P,rho,mu]=ISAAatmosphereModel(h,1) returns T [Rankine, a [ft/s], P
%   [lbf/ft^2], mu [Pa-s] and rho [lbf/ft^3]

% PARAMETERS
% Acceleration due to gravity (m/s^2):
g = 9.80665;
% Ratio of specific heats:
gamma = 1.4;
% Characteristic gas constant (J/Kg/K):
R = 287.0531;
% Lapse rate (K/m):
L = 0.0065;
% Height of troposphere (m):
h_trop = 11000; 
% Height of tropopause (m):
h_strat = 20000; 
% Air density at mean sea level (Kg/m^3):
rho0 = 1.225;
% Ambient pressure at mean sea level (N/m^2):
P0 = 101325;
% Ambient temperature at mean sea level (K):
T0 = 288.15;
% Dynamic viscoity at mean sea level (Pa-s):
mu0=18.27e-6;


if nargin==2 && (nargout==4 || nargout==5) % English units
    hcalc=h*0.3048; %convert to metric
    KtoR=1.8;
    mpsToFps=3.28083989501; % 1 m/s = 
    PaToPsf=0.020885434273; % 1 Pa = # lb/ft^2
    lbDensity=0.06242279606; % 1 kg/m^3 = # lbm/ft^3
elseif nargin==1 && (nargout==4 || nargout==5) % Metric units
    hcalc=h;
    KtoR=1;
    mpsToFps=1;
    PaToPsf=1;
    lbDensity=1;
elseif nargin==2 && nargout==1
    hcalc=h;
    KtoR=1;
    mpsToFps=1;
    PaToPsf=1;
    lbDensity=1;
elseif nargin==3 && nargout==1
    hcalc=h*0.3048; %convert to metric
    KtoR=1.8;
    mpsToFps=3.28083989501; % 1 m/s = 
    PaToPsf=0.020885434273; % 1 Pa = # lb/ft^2
    lbDensity=0.06242279606; % 1 kg/m^3 = # lbm/ft^3
end

if hcalc<=h_trop
    Tcalc=(T0-L*hcalc);
    mu=mu0*(T0+223)/(Tcalc+217)*(Tcalc/(T0))^(1.5)*PaToPsf;
    a=sqrt(gamma*R*Tcalc)*mpsToFps;
    P=(P0*(Tcalc/T0)^((g/(L*R))))*PaToPsf;
    rho=rho0*(Tcalc/T0)^((g/(L*R)-1))*lbDensity;
    T=Tcalc*KtoR;
elseif h_trop<hcalc && hcalc<h_strat
    Tcalc=(T0-L*h_trop);
    mu=mu0*(T0+223)/(Tcalc+217)*(Tcalc/(T0))^(1.5)*PaToPsf;
    a=sqrt(gamma*R*Tcalc)*mpsToFps;
    P=P0*(Tcalc/T0)^(g/(L*R))*exp(-g/R/Tcalc*(hcalc-h_trop))*PaToPsf;
    rho=rho0*(Tcalc/T0)^(g/(L*R))/(Tcalc/T0)*exp(-g/R/Tcalc*(hcalc-h_trop))*lbDensity;
    T=Tcalc*KtoR;
else
    T=T0;
    a=0;
    P=P0;
    rho=rho0;

end

if (nargin==2 || nargin==3) && nargout==1
    M=V/a;
    T=M;
end
end

%Author
%Tyler James Pierce
%tjp644@uw.edu

% Revision History
% Jan/12/2014 - added english units option
% Feb/01/2014 - added dynamic viscosity calculations
%               added mach number calculations
%               


