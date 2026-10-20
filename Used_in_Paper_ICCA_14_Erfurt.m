clear; close all; clc

theta = linspace(0,2*pi,241);
phi   = linspace(0,pi,121);
[theta,phi] = meshgrid(theta,phi);

x = sin(phi).*cos(theta);
y = sin(phi).*sin(theta);
z = cos(phi);

t = linspace(0,1,128)';
blue  = [0.10 0.25 0.75];
white = [1 1 1];
red   = [0.75 0.15 0.10];
cmap = [blue.*(1-t)+white.*t;
        white.*(1-t)+red.*t];

% Four fixed panel positions.
% Titles are added with annotation boxes, not axes titles.
pos = [0.04 0.55 0.39 0.36;
       0.49 0.55 0.39 0.36;
       0.04 0.09 0.39 0.36;
       0.49 0.09 0.39 0.36];

lim = 1.04;


%% Figure 1: four Clifford components, k = 9, n = 5

k = 9;
n = 5;

[F1,F2,F3,F123] = monogenic_components(k,n,theta,phi);

F = {F1,F2,F3,F123};

titles = {sprintf('$[\\widehat{F}_{%d}^{%d}]_{e_1}$',k-1,n), ...
          sprintf('$[\\widehat{F}_{%d}^{%d}]_{e_2}$',k-1,n), ...
          sprintf('$[\\widehat{F}_{%d}^{%d}]_{e_3}$',k-1,n), ...
          sprintf('$[\\widehat{F}_{%d}^{%d}]_{e_{123}}$',k-1,n)};

M = max(abs([F1(:);F2(:);F3(:);F123(:)]));

fig1 = figure('Color','w','Position',[100 100 1200 900]);

for j = 1:4
    ax = axes(fig1,'Position',pos(j,:));

    surf(ax,x,y,z,F{j},'EdgeColor','none');
    daspect(ax,[1 1 1])
    xlim(ax,[-lim lim])
    ylim(ax,[-lim lim])
    zlim(ax,[-lim lim])
    axis(ax,'off')
    view(ax,35,25)
    caxis(ax,[-M M])
    colormap(ax,cmap)

    annotation(fig1,'textbox', ...
        [pos(j,1), pos(j,2)+pos(j,4)+0.012, pos(j,3), 0.04], ...
        'String',titles{j}, ...
        'Interpreter','latex', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',28, ...
        'LineStyle','none');
end

cb = colorbar(ax);
cb.Position = [0.925 0.11 0.018 0.77];
cb.FontSize = 22;
cb.Label.String = 'Clifford coefficient';
cb.Label.FontSize = 24;

exportgraphics(fig1,'Figure1_components_k9_n5.png', ...
    'Resolution',300,'BackgroundColor','white');


%% Figure 2: pointwise norm, k = 9, n = 5

Fnorm = sqrt(F1.^2 + F2.^2 + F3.^2 + F123.^2);

fig2 = figure('Color','w','Position',[100 100 720 650]);
ax = axes(fig2,'Position',[0.06 0.08 0.76 0.80]);

surf(ax,x,y,z,Fnorm,'EdgeColor','none')
daspect(ax,[1 1 1])
xlim(ax,[-lim lim])
ylim(ax,[-lim lim])
zlim(ax,[-lim lim])
axis(ax,'off')
view(ax,35,25)
colormap(ax,turbo)

annotation(fig2,'textbox',[0.06 0.90 0.76 0.05], ...
    'String',sprintf('$|\\widehat{F}_{%d}^{%d}|$',k-1,n), ...
    'Interpreter','latex', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',28, ...
    'LineStyle','none');

cb = colorbar(ax);
cb.Position = [0.87 0.11 0.030 0.76];
cb.FontSize = 22;
cb.Label.String = 'Clifford norm';
cb.Label.FontSize = 24;

exportgraphics(fig2,'Figure2_norm_k9_n5.png', ...
    'Resolution',300,'BackgroundColor','white');


%% Figure 3: fixed k = 9, varying n = 0, 2, 5, 8

k = 9;
nvalues = [0 2 5 8];

norms = cell(1,4);
maxNorm = 0;

for j = 1:4
    n = nvalues(j);
    [A1,A2,A3,A123] = monogenic_components(k,n,theta,phi);
    norms{j} = sqrt(A1.^2 + A2.^2 + A3.^2 + A123.^2);
    maxNorm = max(maxNorm,max(norms{j}(:)));
end

fig3 = figure('Color','w','Position',[100 100 1200 900]);

for j = 1:4
    n = nvalues(j);
    ax = axes(fig3,'Position',pos(j,:));

    surf(ax,x,y,z,norms{j},'EdgeColor','none');
    daspect(ax,[1 1 1])
    xlim(ax,[-lim lim])
    ylim(ax,[-lim lim])
    zlim(ax,[-lim lim])
    axis(ax,'off')
    view(ax,35,25)
    caxis(ax,[0 maxNorm])
    colormap(ax,turbo)

    annotation(fig3,'textbox', ...
        [pos(j,1), pos(j,2)+pos(j,4)+0.012, pos(j,3), 0.04], ...
        'String',sprintf('$|\\widehat{F}_{%d}^{%d}|$',k-1,n), ...
        'Interpreter','latex', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',28, ...
        'LineStyle','none');
end

cb = colorbar(ax);
cb.Position = [0.925 0.11 0.018 0.77];
cb.FontSize = 22;
cb.Label.String = 'Clifford norm';
cb.Label.FontSize = 24;

exportgraphics(fig3,'Figure3_k9_vary_n_norms.png', ...
    'Resolution',300,'BackgroundColor','white');


%% Figure 4: fixed n = 2, varying k = 4, 6, 9, 12

n = 2;
kvalues = [4 6 9 12];

norms = cell(1,4);
maxNorm = 0;

for j = 1:4
    k = kvalues(j);
    [A1,A2,A3,A123] = monogenic_components(k,n,theta,phi);
    norms{j} = sqrt(A1.^2 + A2.^2 + A3.^2 + A123.^2);
    maxNorm = max(maxNorm,max(norms{j}(:)));
end

fig4 = figure('Color','w','Position',[100 100 1200 900]);

for j = 1:4
    k = kvalues(j);
    ax = axes(fig4,'Position',pos(j,:));

    surf(ax,x,y,z,norms{j},'EdgeColor','none');
    daspect(ax,[1 1 1])
    xlim(ax,[-lim lim])
    ylim(ax,[-lim lim])
    zlim(ax,[-lim lim])
    axis(ax,'off')
    view(ax,35,25)
    caxis(ax,[0 maxNorm])
    colormap(ax,turbo)

    annotation(fig4,'textbox', ...
        [pos(j,1), pos(j,2)+pos(j,4)+0.012, pos(j,3), 0.04], ...
        'String',sprintf('$|\\widehat{F}_{%d}^{%d}|$',k-1,n), ...
        'Interpreter','latex', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',28, ...
        'LineStyle','none');
end

cb = colorbar(ax);
cb.Position = [0.925 0.11 0.018 0.77];
cb.FontSize = 22;
cb.Label.String = 'Clifford norm';
cb.Label.FontSize = 24;

exportgraphics(fig4,'Figure4_vary_k_n2_norms.png', ...
    'Resolution',300,'BackgroundColor','white');


function [F1,F2,F3,F123] = monogenic_components(k,n,theta,phi)

    if n == 0
        Pk  = jacobi_poly(k,0,0,cos(phi));
        dPk = jacobi_poly_derivative(k,0,0,cos(phi));

        c0 = 1/sqrt(8*pi);

        F1 = c0*sin(phi).*(k*Pk - cos(phi).*dPk).*cos(theta);
        F2 = c0*sin(phi).*(k*Pk - cos(phi).*dPk).*sin(theta);
        F3 = c0*(k*cos(phi).*Pk + sin(phi).^2.*dPk);
        F123 = zeros(size(phi));

        normF = sqrt(k/2);

    else
        Pkn  = jacobi_poly(k-n,   n,n,cos(phi));
        Pkn1 = jacobi_poly(k-n-1, n,n,cos(phi));

        pref12 = sin(phi).^(n-1)/(2^(n+1)*sqrt(pi));
        bracket = (k-n)*Pkn - k*cos(phi).*Pkn1;

        F1 = pref12.*bracket.*cos((n+1)*theta);
        F2 = pref12.*bracket.*sin((n+1)*theta);

        pref3 = k*sin(phi).^n/(2^(n+1)*sqrt(pi));
        F3    =  pref3.*Pkn1.*cos(n*theta);
        F123  = -pref3.*Pkn1.*sin(n*theta);

        normF = exp(0.5*(2*gammaln(k+1) ...
                    - gammaln(k-n) ...
                    - gammaln(k+n+1)));
    end

    F1   = F1/normF;
    F2   = F2/normF;
    F3   = F3/normF;
    F123 = F123/normF;
end


function P = jacobi_poly(N,alpha,beta,x)

    if N < 0
        P = zeros(size(x));
        return
    end

    P = zeros(size(x));

    for s = 0:N
        c = nchoosek(N+alpha,N-s)*nchoosek(N+beta,s)/2^N;
        P = P + c.*(x-1).^s.*(x+1).^(N-s);
    end
end


function dP = jacobi_poly_derivative(N,alpha,beta,x)

    if N <= 0
        dP = zeros(size(x));
        return
    end

    dP = zeros(size(x));

    for s = 0:N
        c = nchoosek(N+alpha,N-s)*nchoosek(N+beta,s)/2^N;

        if s >= 1
            dP = dP + c*s.*(x-1).^(s-1).*(x+1).^(N-s);
        end

        if N-s >= 1
            dP = dP + c*(N-s).*(x-1).^s.*(x+1).^(N-s-1);
        end
    end
end
