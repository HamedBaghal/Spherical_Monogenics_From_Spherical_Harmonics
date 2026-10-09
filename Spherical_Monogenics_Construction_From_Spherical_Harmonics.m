clear; close all; clc

k = input('Enter homogeniety degree k to produce all spgerical monogenics from 0 to k-1: ');

theta = linspace(0,2*pi,241);
phi   = linspace(0,pi,121);
[theta,phi] = meshgrid(theta,phi);

x = sin(phi).*cos(theta);
y = sin(phi).*sin(theta);
z = cos(phi);

% blue-white-red colormap
t = linspace(0,1,128)';
blue  = [0.10 0.25 0.75];
white = [1 1 1];
red   = [0.75 0.15 0.10];
cmap = [blue.*(1-t)+white.*t; white.*(1-t)+red.*t];

for n = 0:k-1

    if n == 0
        Pk  = jacobiP(k,0,0,cos(phi));
        dPk = (k+1)/2 * jacobiP(k-1,1,1,cos(phi));

        F1   = sin(phi).*(k*Pk-cos(phi).*dPk).*cos(theta)/sqrt(8*pi);
        F2   = sin(phi).*(k*Pk-cos(phi).*dPk).*sin(theta)/sqrt(8*pi);
        F3   = (k*cos(phi).*Pk+sin(phi).^2.*dPk)/sqrt(8*pi);
        F123 = zeros(size(phi));

        normF = sqrt(k/2);

    else
        Pkn  = jacobiP(k-n,n,n,cos(phi));
        Pkn1 = jacobiP(k-n-1,n,n,cos(phi));

        A = sin(phi).^(n-1)/(2^(n+1)*sqrt(pi));
        B = (k-n)*Pkn-k*cos(phi).*Pkn1;

        F1 = A.*B.*cos((n+1)*theta);
        F2 = A.*B.*sin((n+1)*theta);

        C = k*sin(phi).^n/(2^(n+1)*sqrt(pi));
        F3   =  C.*Pkn1.*cos(n*theta);
        F123 = -C.*Pkn1.*sin(n*theta);

        normF = factorial(k)/sqrt(factorial(k-n-1)*factorial(k+n));
    end

    F1   = F1/normF;
    F2   = F2/normF;
    F3   = F3/normF;
    F123 = F123/normF;

    F = {F1,F2,F3,F123};
    names = {'e_1','e_2','e_3','e_{123}'};
    files = {'e1','e2','e3','e123'};

    M = max(abs([F1(:);F2(:);F3(:);F123(:)]));

    % four separate component plots
    for j = 1:4
        figure('Color','w','Position',[100 100 560 470])
        surf(x,y,z,F{j},'EdgeColor','none')
        axis equal off
        view(35,25)
        caxis([-M M])
        colormap(cmap)
        colorbar
        title(sprintf('$[\\widehat{F}_{%d}^{%d}]_{%s}$',k-1,n,names{j}), ...
              'Interpreter','latex','FontSize',20)
        exportgraphics(gcf,sprintf('Fhat_k%d_n%d_%s.png',k,n,files{j}), ...
                       'Resolution',300)
    end

    % all four components together
    figure('Color','w','Position',[100 100 1000 820])
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact')
    for j = 1:4
        ax = nexttile;
        surf(ax,x,y,z,F{j},'EdgeColor','none')
        axis(ax,'equal','off')
        view(ax,35,25)
        caxis(ax,[-M M])
        colormap(ax,cmap)
        title(ax,sprintf('$[\\widehat{F}_{%d}^{%d}]_{%s}$',k-1,n,names{j}), ...
              'Interpreter','latex','FontSize',18)
    end
    cb = colorbar;
    cb.Layout.Tile = 'east';
    exportgraphics(gcf,sprintf('Fhat_k%d_n%d_components.png',k,n), ...
                   'Resolution',300)

    % pointwise Clifford norm
    Fnorm = sqrt(F1.^2+F2.^2+F3.^2+F123.^2);

    figure('Color','w','Position',[100 100 560 470])
    surf(x,y,z,Fnorm,'EdgeColor','none')
    axis equal off
    view(35,25)
    colormap(turbo)
    colorbar
    title(sprintf('$|\\widehat{F}_{%d}^{%d}|$',k-1,n), ...
          'Interpreter','latex','FontSize',20)
    exportgraphics(gcf,sprintf('Fhat_k%d_n%d_norm.png',k,n), ...
                   'Resolution',300)
end
