clear; close all; clc

theta = linspace(0,2*pi,181);
phi   = linspace(0,pi,91);
[theta,phi] = meshgrid(theta,phi);

x = sin(phi).*cos(theta);
y = sin(phi).*sin(theta);
z = cos(phi);

a = sqrt(5)/(8*sqrt(pi));
b = 3*sqrt(5)/(4*sqrt(pi));

F1   = -a*sin(phi).^3.*cos(3*theta);
F2   = -a*sin(phi).^3.*sin(3*theta);
F3   =  b*sin(phi).^2.*cos(phi).*cos(2*theta);
F123 = -b*sin(phi).^2.*cos(phi).*sin(2*theta);

F = {F1,F2,F3,F123};
titles = {'$[\widehat F_3^2]_{e_1}$', ...
          '$[\widehat F_3^2]_{e_2}$', ...
          '$[\widehat F_3^2]_{e_3}$', ...
          '$[\widehat F_3^2]_{e_{123}}$'};

M = max(abs([F1(:);F2(:);F3(:);F123(:)]));

t = linspace(0,1,128)';
blue  = [0.10 0.25 0.75];
white = [1 1 1];
red   = [0.75 0.15 0.10];

cmap = [blue.*(1-t)+white.*t;
        white.*(1-t)+red.*t];

figure('Color','w','Position',[100 100 950 760])
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for j = 1:4
    ax = nexttile;
    surf(ax,x,y,z,F{j},'EdgeColor','none');
    axis(ax,'equal','off')
    view(ax,35,25)
    caxis(ax,[-M M])
    colormap(ax,cmap)
    title(ax,titles{j},'Interpreter','latex','FontSize',20)
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Clifford coefficient';

exportgraphics(gcf,'Figure2_components.png','Resolution',300);
