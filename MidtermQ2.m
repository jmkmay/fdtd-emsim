clc;
clear all;

epsilon_0 = 8.854e-12 ;
mu_0 = (4*pi)*10^-7;


A = 1;

[R T] = meshgrid(0:1:1000,0:pi/30:2*pi);

% 
% theta = linspace(0,2*pi);
% r = linspace(0,1000);

S = sqrt(epsilon_0/mu_0)*A^2*(sin(T)).^2/(2*R.^2);

x = R.cos(T);
y = R.sin(T);

h = polar(x,y);
counterf(x,y,S);