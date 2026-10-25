%%%% Model obiektu rozpisany z użyciem solvera ode45
clear; clc; close all;

%% Parametry
p.A2 = 430; p.C1 = 0.85; p.a1 = 18; p.a2 = 24; p.tau = 70;
F1_0 = 55; FD_0 = 10;

%% Punkt pracy
h1_0 = ((F1_0 + FD_0)/p.a1)^2;
h2_0 = ((F1_0 + FD_0)/p.a2)^2;

%% Wejścia (te same zmiany co w scenariuszu DMC): skoki F1in i FD o +/-50%
F1in = @(t) F1_0*((t < 300) + 1.5*(t >= 300 & t < 1200) + 0.5*(t >= 1200));
FD   = @(t) FD_0*((t < 700) + 0.5*(t >= 700 & t < 2000) + 1.5*(t >= 2000));
F1   = @(t) F1in(t - p.tau);

%% Rozwiązanie
opts = odeset('RelTol',1e-6, 'AbsTol',1e-8, 'MaxStep',1);
[t, h] = ode45(@(t,h) zbiorniki(t,h,p,F1,FD), [0 3000], [h1_0; h2_0], opts);

%% Wykresy
folder = fullfile(fileparts(mfilename('fullpath')), 'zad1_plots');
f = figure('Color', 'w', 'Visible', 'off');
subplot(2,1,1); plot(t, h(:,2)); grid on; ylabel('h_2 [cm]');
subplot(2,1,2); plot(t, F1in(t)); hold on; plot(t, FD(t), '--'); hold off;
grid on; ylabel('F [cm^3/s]'); xlabel('t [s]');
legend('F_{1in}', 'F_D', 'Location', 'best');
exportgraphics(f, fullfile(folder, 'zad1_model_ode45.png'), 'Resolution', 300); close(f);

%% Model
function dh = zbiorniki(t, h, p, F1, FD)
h1 = max(h(1), 1e-6);
h2 = max(h(2), 0);
F2 = p.a1*sqrt(h1);
F3 = p.a2*sqrt(h2);
dh = [ (F1(t) + FD(t) - F2) / (2*p.C1*h1) ;
    (F2 - F3) / p.A2 ];
end