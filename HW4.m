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

%% Initial guesses for parameters
param.kf = -0.000000000000000000000000000032; % Adjust this value
kf = param.kf;

param.ki = 0.001; % Adjust this value
ki = param.ki;

param.alpha = 0.0005; % Adjust this value
alpha = param.alpha;

param.beta = -0.00001; % Adjust this value
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

% Adjust kf, ki, alpha, and beta to meet requirements
% Use the stepinfo function to get overshoot and other metrics
% Use the bandwidth function to get the bandwidth of the system


info = stepinfo(G_close_loop);
overshoot = info.Overshoot;
bw = bandwidth(G_close_loop);

fprintf('Overshoot: %f\n', overshoot);
fprintf('Bandwidth: %f rad/sec\n', bw);

% Iterate to adjust parameters based on the results
% Example adjustments (these would typically be done manually based on the output)
% param.kf = new_kf;
% param.ki = new_ki;
% param.alpha = new_alpha;
% param.beta = new_beta;
