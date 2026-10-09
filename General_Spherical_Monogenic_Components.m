clear; close all; clc

% Three-dimensional case only: R_3.
% k = harmonic homogeneity degree, so the monogenic degree is k-1.
k = input('Harmonic homogeneity degree k (k >= 1): ');
n = input('Basis index n, with 0 <= n <= k-1: ');

if k < 1 || k ~= floor(k)
    error('k must be an integer with k >= 1.');
end
if n < 0 || n > k-1 || n ~= floor(n)
    error('n must be an integer with 0 <= n <= k-1.');
end

% true  -> plots the normalized basis element \widehat{F}_{k-1}^n
% false -> plots the raw basis element F_{k-1}^n
normalize = true;

outDir = sprintf('monogenic_components_k%d_n%d',k,n);
if ~exist(outDir,'dir')
    mkdir(outDir)
end

theta = linspace(0,2*pi,241);
phi   = linspace(0,pi,121);
[theta,phi] = meshgrid(theta,phi);

x = sin(phi).*cos(theta);
y = sin(phi).*sin(theta);
z = cos(phi);

% ---------------------------------------------------------------
% Raw Clifford coefficients.
% For n >= 1 these are the unsimplified formulas (3.35)--(3.38).
% For n = 0 they are obtained directly from equation (3.12).
% ---------------------------------------------------------------

if n == 0
    Pk  = jacobi_poly(k,0,0,cos(phi));
    dPk = jacobi_poly_derivative(k,0,0,cos(phi));

    c0 = 1/sqrt(8*pi);

    F1 = c0*sin(phi) .* ...
        (k*Pk - cos(phi).*dPk) .* cos(theta);

    F2 = c0*sin(phi) .* ...
        (k*Pk - cos(phi).*dPk) .* sin(theta);

    F3 = c0*(k*cos(phi).*Pk + sin(phi).^2.*dPk);

    F123 = zeros(size(phi));

else
    Pkn  = jacobi_poly(k-n,   n, n, cos(phi));
    Pkn1 = jacobi_poly(k-n-1, n, n, cos(phi));

    pref12 = sin(phi).^(n-1)/(2^(n+1)*sqrt(pi));
    bracket = (k-n)*Pkn - k*cos(phi).*Pkn1;

    F1 = pref12 .* bracket .* cos((n+1)*theta);
    F2 = pref12 .* bracket .* sin((n+1)*theta);

    pref3 = k*sin(phi).^n/(2^(n+1)*sqrt(pi));
    F3    =  pref3 .* Pkn1 .* cos(n*theta);
    F123  = -pref3 .* Pkn1 .* sin(n*theta);
end

% ---------------------------------------------------------------
% Normalize using Proposition 3.6.
% ---------------------------------------------------------------

if normalize
    if n == 0
        normF = sqrt(k/2);
    else
        logNormSq = 2*gammaln(k+1) ...
                  - gammaln(k-n) ...
                  - gammaln(k+n+1);
        normF = exp(0.5*logNormSq);
    end

    F1   = F1/normF;
    F2   = F2/normF;
    F3   = F3/normF;
    F123 = F123/normF;

    fileTag = 'normalized';
else
    fileTag = 'raw';
end

F = {F1,F2,F3,F123};
componentNames = {'e_1','e_2','e_3','e_{123}'};
fileNames = {'e1','e2','e3','e123'};

% LaTeX titles.  Braces around F make MATLAB display the widehat correctly.
if normalize
    componentTitles = cell(1,4);
    for j = 1:4
        componentTitles{j} = sprintf( ...
            '$[\\widehat{F}_{%d}^{%d}]_{%s}$', ...
            k-1,n,componentNames{j});
    end
    normTitle = sprintf('$|\\widehat{F}_{%d}^{%d}|$',k-1,n);
else
    componentTitles = cell(1,4);
    for j = 1:4
        componentTitles{j} = sprintf( ...
            '$[F_{%d}^{%d}]_{%s}$', ...
            k-1,n,componentNames{j});
    end
    normTitle = sprintf('$|F_{%d}^{%d}|$',k-1,n);
end

M = max(abs([F1(:);F2(:);F3(:);F123(:)]));

% Blue-white-red map for signed Clifford coefficients.
t = linspace(0,1,128)';
blue  = [0.10 0.25 0.75];
white = [1 1 1];
red   = [0.75 0.15 0.10];
cmap = [blue.*(1-t)+white.*t;
        white.*(1-t)+red.*t];

% ---------------------------------------------------------------
% Four separate component figures.
% ---------------------------------------------------------------

% for j = 1:4
%     fig = figure('Color','w','Position',[100 100 560 470]);
%     surf(x,y,z,F{j},'EdgeColor','none');
%     axis equal off
%     view(35,25)
%     caxis([-M M])
%     colormap(cmap)
%     cb = colorbar;
%     cb.Label.String = 'Clifford coefficient';
%     title(componentTitles{j},'Interpreter','latex','FontSize',20)
% 
%     exportgraphics(fig,fullfile(outDir, ...
%         sprintf('%s_component_%s_k%d_n%d.png', ...
%         fileTag,fileNames{j},k,n)), ...
%         'Resolution',300);
% end

% ---------------------------------------------------------------
% All four components in one figure.
% ---------------------------------------------------------------

figAll = figure('Color','w','Position',[100 100 1000 820]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for j = 1:4
    ax = nexttile;
    surf(ax,x,y,z,F{j},'EdgeColor','none');
    axis(ax,'equal','off')
    view(ax,35,25)
    caxis(ax,[-M M])
    colormap(ax,cmap)
    title(ax,componentTitles{j},'Interpreter','latex','FontSize',18)
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Clifford coefficient';

exportgraphics(figAll,fullfile(outDir, ...
    sprintf('%s_all_components_k%d_n%d.png',fileTag,k,n)), ...
    'Resolution',300);

% ---------------------------------------------------------------
% Pointwise Clifford absolute value, computed from the components.
% No simplified closed formula is used here.
% ---------------------------------------------------------------

Fnorm = sqrt(F1.^2 + F2.^2 + F3.^2 + F123.^2);

figNorm = figure('Color','w','Position',[100 100 560 470]);
surf(x,y,z,Fnorm,'EdgeColor','none');
axis equal off
view(35,25)
colormap(turbo)
cb = colorbar;
cb.Label.String = 'Clifford norm';
title(normTitle,'Interpreter','latex','FontSize',20)

exportgraphics(figNorm,fullfile(outDir, ...
    sprintf('%s_absolute_value_k%d_n%d.png',fileTag,k,n)), ...
    'Resolution',300);

fprintf('\nSaved figures in: %s\n',outDir);
fprintf('  4 separate component figures\n');
fprintf('  1 combined 2-by-2 component figure\n');
fprintf('  1 pointwise absolute-value figure\n');


function P = jacobi_poly(N,alpha,beta,x)
% Jacobi polynomial P_N^{(alpha,beta)}(x), using the normalization
% stated in the paper. N, alpha and beta are nonnegative integers.

    if N < 0
        P = zeros(size(x));
        return
    end

    P = zeros(size(x));
    for s = 0:N
        c = nchoosek(N+alpha,N-s) * nchoosek(N+beta,s) / 2^N;
        P = P + c .* (x-1).^s .* (x+1).^(N-s);
    end
end


function dP = jacobi_poly_derivative(N,alpha,beta,x)
% Derivative with respect to the scalar Jacobi argument.
% This differentiates the explicit Jacobi sum directly.

    if N <= 0
        dP = zeros(size(x));
        return
    end

    dP = zeros(size(x));

    for s = 0:N
        c = nchoosek(N+alpha,N-s) * nchoosek(N+beta,s) / 2^N;

        if s >= 1
            dP = dP + c*s .* ...
                (x-1).^(s-1) .* (x+1).^(N-s);
        end

        if N-s >= 1
            dP = dP + c*(N-s) .* ...
                (x-1).^s .* (x+1).^(N-s-1);
        end
    end
end
