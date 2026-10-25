%%%% Model obiektu rozpisany z użyciem solvera ode45 (równania: zad1_obiekt.ode)

clear; clc; close all;

%% Parametry i punkt pracy
[model, F1_0, FD_0, h1_0, h2_0] = zad1_obiekt.parametry();
folder = fullfile(fileparts(mfilename('fullpath')), 'zad1_plots');

%% Wejścia
F1in = @(tt) zad1_scenariusze.przebieg_F1(tt, F1_0);
FD = @(tt) zad1_scenariusze.przebieg_FD(tt, FD_0);
F1 = @(tt) F1in(tt - model.tau);

%% Rozwiązanie
opts = odeset('RelTol',1e-6, 'AbsTol',1e-8, 'MaxStep',1);
[t, h] = ode45(@(tt,h) zad1_obiekt.ode(tt, h, model, F1, FD), [0 3000], [h1_0; h2_0], opts);

%% Wykresy
f = zad1_wykresy.nowy();
subplot(2,1,1); plot(t, h(:,2)); grid on; ylabel('h_2 [cm]');
subplot(2,1,2); plot(t, F1in(t)); hold on; plot(t, FD(t), 'black'); hold off;
grid on; ylabel('F [cm^3/s]'); xlabel('t [s]');
legend('F_{1in}', 'F_D', 'Location', 'best');
zad1_wykresy.zapisz(f, fullfile(folder, 'zad1_model_ode45.png'));
