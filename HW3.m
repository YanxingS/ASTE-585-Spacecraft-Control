clc; clear; close all;

%% HW 3 problem 1
% problem wants us to design a PD controller for a momentum biased
% (negative y axis) satellite, to do this we first solve for kpx and kpz

%% Define parameters 
global param
param.w0 = 7.326e-5;
w0 = param.w0;

param.Hb = 30;
Hb = param.Hb;

param.Tdx = 6e-5;
param.Tdz = 4e-5;
Tdx = param.Tdx;
Tdz = param.Tdz;

param.phi_ss = deg2rad(0.4);
param.psi_ss = deg2rad(0.4);
phi_ss = param.phi_ss;
psi_ss = param.psi_ss;

param.Ix = 90;
param.Iy = 100;
param.Iz = 110;
Ix = param.Ix;
Iy = param.Iy;
Iz = param.Iz;

param.damping = 0.7; % an assumption made for damping coefficient
damping = param.damping;
%% First we solve for kpx and kpz

syms kpx kpz

eqn1 = (param.phi_ss/param.Tdx) == (kpz+param.w0*param.Hb)/(kpz*kpx+(param.w0*param.Hb)^2+param.w0*param.Hb*(kpx+kpz));
eqn2 = (param.psi_ss/param.Tdz) == (kpx+param.w0*param.Hb)/(kpz*kpx+(param.w0*param.Hb)^2+param.w0*param.Hb*(kpx+kpz));

solution  = solve([eqn1,eqn2],[kpx kpz]);

kpx_sol = vpa(solution.kpx);
kpz_sol = vpa(solution.kpz);

clear kpx kpz 

param.kpx = kpx_sol(1);
param.kpz = kpz_sol(2);
kpx = param.kpx;
kpz = param.kpz;

%% Solving Systems of equation to determine kdz kdx wn1 wn2

fun = @Momentum_Biased_Satellite;
x0 = [22 12 0 0];
x = fsolve(fun,x0);

param.kdx = x(1);
param.kdz = x(2);
param.wn1 = x(3);
param.wn2 = x(4);
kdx = param.kdx;
kdz = param.kdz;
wn1 = param.wn1;
wn2 = param.wn2;
%% Finding TF and create bode plot 

% phi/Tdx
numerator_phi_Tdx = [1, 5646000545634665/18014398509481984, -0.0000000000000000000036943075969773469081659352796609];
denominator = [1, 1746502305688195/72057594037927936, 0.00028550509825247877432026030549501, -0.00000276145032430467798943844768103, -0.0000000000000000000003492525853243666840813208344817];
phi_Tdx_TF = tf(numerator_phi_Tdx,denominator);
% phi/Tdz
numerator_phi_Tdz = [0.3333,0];
phi_Tdz_TF = tf(numerator_phi_Tdz,denominator);
% psi/Tdx
numerator_psi_Tdx = [0.2727,0];
psi_Tdx_TF = tf(numerator_psi_Tdx,denominator);
% psi/Tdz
numerator_psi_Tdz = [1, 5646000545634665/18014398509481984, -0.0000000000000000000036943075969773469081659352796609];
psi_Tdz_TF = tf(numerator_psi_Tdz,denominator);

% create bode plot

bode(phi_Tdx_TF,phi_Tdz_TF,psi_Tdx_TF,psi_Tdz_TF);
legend('phi_Tdx','phi_Tdz','psi_Tdx','psi_Tdz');

%% problem 1b satellite now only has roll measurement
% first we need to find the necessary momentum bia wheel Hb
clear Hb kpx damping

% Initial guess for a, damping
param.damping = 0.7;
damping = param.damping;

param.a = 0.5;
a = param.a;

param.Hb = (param.Tdz-param.a*Tdx)/(param.w0*param.psi_ss);
param.kpx = (param.Tdx/param.psi_ss)-(param.w0*param.Hb);
Hb = param.Hb;
kpx = param.kpx;
%% Solving for Kdx

clear fun x0 x

fun = @roll_only_momentum_bias;
x0 = [22 0 0];
x = fsolve(fun,x0);
param.kdx = x(1);
param.wn1 = x(2);
param.wn2 = x(3);
kdx = param.kdx;
wn1 = param.wn1;
wn2 = param.wn2;

%% Finding Transfer function and create bode plot

denominator = [1,kdx*Iz...
    ,w0*Hb*(Iz+Ix)+Hb^2+(kpx*Iz+a*Hb*kdx)...
    ,a*Hb*kpx+kdx*w0*Hb...
    ,w0*Hb*(w0*Hb+kpx)];

%phi/Tdx
numerator_phi_Tdx = [1,0,1.3022e-05];
phi_Tdx_TF = tf(numerator_phi_Tdx,denominator);
%phi/Tdz
numerator_phi_Tdz = [-Hb/Ix,0];
phi_Tdz_TF = tf(numerator_phi_Tdz,denominator);
%psi/Tdx
numerator_psi_Tdx = [(Hb+a*kdx)/Iz,a*kpx/Iz];
psi_Tdx_TF = tf(numerator_psi_Tdx,denominator);
%psi/Tdz
numerator_psi_Tdz = [1,kdx/Ix,(kpx+w0*Hb)/Ix];
psi_Tdz_TF = tf(numerator_psi_Tdz,denominator);

% create bode plot
figure(2)
bode(phi_Tdx_TF,phi_Tdz_TF,psi_Tdx_TF,psi_Tdz_TF);
legend('phi_Tdx_TF','phi_Tdz_TF','psi_Tdx_TF','psi_Tdz_TF');
%% Function that finds kdx, kdz, wn1, wn2

function F = Momentum_Biased_Satellite(x)
% definition of x = [kdx kdz wn1 wn2]
global param

Ix = param.Ix;
Iy = param.Iy;
Iz = param.Iz;

kpx = param.kpx;
kpz = param.kpz;

damping = param.damping;

Hb = param.Hb;
w0 = param.w0;

F(1) = 2*(damping*x(3)+damping*x(4))-(x(1)*Iz+x(2)*Ix)/(Ix*Iz);
F(2) = x(3)^2+x(4)^2+4*damping*damping*x(3)*x(4)-((x(1)*x(2)+Hb^2+(kpx*Iz+w0*Hb*Iz+kpz*Ix+w0*Hb*Ix))/(Ix*Iz));
F(3) = 2*x(3)*x(4)*(damping*x(4)+damping*x(3))-(x(2)*(kpz+w0*Hb)+(x(1)*(kpx+w0*Hb)))/(Ix*Iz);
F(4) = x(3)^2+x(4)^2-((kpx+w0*Hb)*(kpz+w0*Hb))/(Ix*Iz);
end

%% Function that computes momentum wheel bias for 1b 

function F = roll_only_momentum_bias(x)
% x = [kdx wn1 wn2]
global param

Ix = param.Ix;
Iy = param.Iy;
Iz = param.Iz;

kpx = param.kpx;
kpz = param.kpz;

damping = param.damping;

Hb = param.Hb;
w0 = param.w0;
a = param.a;

F(1) = 2*(damping*x(2)+damping*x(3))-(x(1)*Iz)/(Ix*Iz);
F(2) = x(2)^2+x(3)^2+4*damping*damping*x(2)*x(3)-(w0*Hb*(Ix+Iz)+Hb^2+(kpx*Iz+a*Hb*x(1)))/(Ix*Iz);
F(3) = 2*x(2)*x(3)*(damping*x(3)+damping*x(2));
end