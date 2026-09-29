clear; close all; clc

% Figure 1 in the paper: k = 4, n = 2

theta = linspace(0,2*pi,181);
phi   = linspace(0,pi,91);
[theta,phi] = meshgrid(theta,phi);

x = sin(phi).*cos(theta);
y = sin(phi).*sin(theta);
z = cos(phi);

% H_4^{cos,2} and H_4^{sin,2}
% P_2^{(2,2)}(t) = 7t^2 - 1

A = sin(phi).^2 .* (7*cos(phi).^2 - 1) /(8*sqrt(pi));

Hcos = A.*cos(2*theta);
Hsin = A.*sin(2*theta);

% |F_3^2| from the explicit formula for the raw monogenic

Fnorm = sin(phi).^2 .* sqrt(1 + 35*cos(phi).^2) /(4*sqrt(pi));

% Blue-white-red colour map for the signed harmonic modes

t = linspace(0,1,128)';
blue  = [0.10 0.25 0.75];
white = [1 1 1];
red   = [0.75 0.15 0.10];

cmap = [blue.*(1-t) + white.*t;
        white.*(1-t) + red.*t];

M = max(abs([Hcos(:); Hsin(:)]));

% H_4^{cos,2}

figure('Color','w','Position',[100 100 500 420])
surf(x,y,z,Hcos,'EdgeColor','none')
axis equal off
view(35,25)
caxis([-M M])
colormap(cmap)
colorbar
exportgraphics(gcf,'construction_cos_k4_n2.png','Resolution',300)

% H_4^{sin,2}

figure('Color','w','Position',[100 100 500 420])
surf(x,y,z,Hsin,'EdgeColor','none')
axis equal off
view(35,25)
caxis([-M M])
colormap(cmap)
colorbar
exportgraphics(gcf,'construction_sin_k4_n2.png','Resolution',300)

% |F_3^2|

figure('Color','w','Position',[100 100 500 420])
surf(x,y,z,Fnorm,'EdgeColor','none')
axis equal off
view(35,25)
colormap(turbo)
colorbar
exportgraphics(gcf,'construction_norm_k4_n2.png','Resolution',300)
