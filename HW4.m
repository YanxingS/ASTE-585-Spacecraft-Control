clc; clear; close all;

%% HW 4 ASTE 585

% design an attitude controller for the pitch axis of a satellite 

%% Parameter definition

global param

param.Ix = 90;
param.Iy = 80;
param.Iz = 100;

Ix = param.Ix;
Iy = param.Iy;
Iz = param.Iz; 

param.Hb = 27;
Hb = param.Hb;

param.w0 = 2*pi/(60*60*24);
w0 = param.w0;

param.epsilon = w0*(1-1/sqrt(1+(w0*Iz/Hb)));
epsilon = param.epsilon;

param.wn = sqrt(Hb^2/(Ix*Iz));
wn = param.wn;

%% Initial guesses for parameters
param.kf = 0.00000006; % Adjust this value
kf = param.kf;
    
param.ki = 0.013; % Adjust this value
ki = param.ki;

param.alpha = 0.56; % Adjust this value
alpha = param.alpha;

param.beta = 0.6; % Adjust this value
beta = param.beta;

s = tf('s');
K_s = -kf*(1+(ki/s))*((s+alpha)/(s+beta));

%% define pitch transfer function (theta/Tcy)

Ix = param.Ix;
Iz = param.Iz;
w0 = param.w0;

G_theta_s_num = 1;
G_theta_s_den = [Iy,0,0];
G_theta_s = tf(G_theta_s_num,G_theta_s_den);

%% define close loop transfer function 

H = 1; % perfect feedback
G_open_loop = K_s*G_theta_s;
G_close_loop = feedback(G_open_loop,H);

% Step response analysis
figure(1);
step(G_close_loop);
stepinfo(G_close_loop)

% Bode plot analysis
figure(2);
bode(G_close_loop);

% root locus plot analysis
figure(3);
rlocus(G_open_loop);
% Adjust kf, ki, alpha, and beta to meet requirements
% Use the stepinfo function to get overshoot and other metrics
% Use the bandwidth function to get the bandwidth of the system


info = stepinfo(G_close_loop);
overshoot = info.Overshoot;
bw = bandwidth(G_close_loop);

%% define transfer function for roll phi/gamma anad yaw psi/gamma
determinant = Ix*Iz*(s^2+w0^2)*(s^2+wn^2);
G_phi_s = (-Hb^2*(s^2+(w0-epsilon)^2))/determinant; %phi/gamma
G_psi_s = (s*Hb*(1+G_phi_s))/(s^2+w0*Hb);

%% define controller for roll and pitch 

K_s_roll = -((1+0.013/s)*(s+0.05)/(s+0.65));

%% define close loop transfer function for roll and yaw
opts = stepDataOptions('StepAmplitude',deg2rad(5)); % adjust step input magnitude

G_cl_phi = feedback(K_s_roll*G_phi_s,H);
figure(4)
step(G_cl_phi,opts);

G_cl_psi = feedback(K_s_roll*G_psi_s,H);
figure(5)
step(G_cl_psi,opts);

%% define transfer function for disturbance input 
G_x = (Iz*s^2+w0*Hb)/(determinant);
G_z = (Hb*s)/(determinant);
    
%% define disturbance equation and perform laplace transformation
syms t 

T_di_t = 10e-5+10e-4*sin(w0*t);

T_di_x = laplace(T_di_t);

%% define actuator parameters and TF

param.tao_gamma = 5;
tao_gamma = param.tao_gamma;

param.tao_w = 5;
tao_w = param.tao_w;

roll_motor = 1/(1+tao_gamma*s);
pitch_motor = 1/(1+tao_w*s);

%% construct new close loop with actuator model
G_cl_phi = feedback(K_s_roll*roll_motor*G_phi_s,H);
G_cl_theta = feedback(K_s*G_theta_s*pitch_motor,H);

figure(6)
bode(G_cl_phi);
figure(7)
bode(G_cl_theta);


