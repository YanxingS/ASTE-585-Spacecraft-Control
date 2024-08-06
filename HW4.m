clc; clear; close all;

%% HW 4 ASTE 585

% design an attitude controller for the pitch axis of a satellite 

%% Parameter definition

global param

param.Ix = 90;
param.Iy = 80;
param.Ix = 100;

Ix = param.Ix;
Iy = param.Iy;
Iz = param.Iz; 

param.Hb = 27;
Hb = param.Hb;

param.w0 = 2*pi/(60*60*24);
w0 = param.w0;

%% define controller transfer function and parameters

param.kf = 0.5;
kf = param.kf;

param.ki = 0.01;
ki = param.ki;

param.alpha = 0.004;
alpha = param.alpha;

param.beta = 0.0005;
beta = param.beta;

s = tf('s');
K_s = -kf*(1+ki/s)*((s+alpha)/(s+beta));

% K_s_num = [-kf,-(kf*alpha+kf*ki),-ki*kf*alpha];
% K_s_den = [1,beta,0];
% K_s_1 = tf(K_s_num,K_s_den);

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
3
step(G_close_loop);