%%%% Model obiektu rozpisany metodą Eulera
clear; clc; close all;

%% Parametry
p.A2 = 430; p.C1 = 0.85; p.a1 = 18; p.a2 = 24; p.tau = 70;
F1_0 = 55; FD_0 = 10;
h1_0 = ((F1_0 + FD_0)/p.a1)^2;
h2_0 = ((F1_0 + FD_0)/p.a2)^2;

%% Dyskretyzacja
Ts = 0.5;                 % krok [s]
Tk = 3000;                % czas symulacji [s]
N  = round(Tk/Ts);
d  = round(p.tau/Ts);     % opóźnienie w próbkach

t    = (0:N)*Ts;
h1   = zeros(1,N+1);  h2 = zeros(1,N+1);
h1(1) = h1_0;  h2(1) = h2_0;

%% Wejścia (te same zmiany co w scenariuszu DMC): skoki F1in i FD o +/-50%
F1in = F1_0*ones(1,N+1);
F1in(t >= 300)  = 1.5*F1_0;
F1in(t >= 1200) = 0.5*F1_0;
FD = FD_0*ones(1,N+1);
FD(t >= 700)  = 0.5*FD_0;
FD(t >= 2000) = 1.5*FD_0;

%% Pętla symulacji
for k = 1:N
    if k > d
        F1 = F1in(k-d);
    else
        F1 = F1_0;
    end

    F2 = p.a1*sqrt(h1(k));
    F3 = p.a2*sqrt(h2(k));

    % Równanie różnicowe

    h1(k+1) = h1(k) + Ts*(F1 + FD(k) - F2)/(2*p.C1*h1(k));
    h2(k+1) = h2(k) + Ts*(F2 - F3)/p.A2;

    h1(k+1) = max(h1(k+1), 1e-6);
    h2(k+1) = max(h2(k+1), 0);
end

%% Wykresy
folder = fullfile(fileparts(mfilename('fullpath')), 'zad1_plots');
f = figure('Color', 'w', 'Visible', 'off');
subplot(2,1,1); plot(t, h2); grid on; ylabel('h_2 [cm]');
subplot(2,1,2); stairs(t, F1in); hold on; stairs(t, FD, '--'); hold off;
grid on; ylabel('F [cm^3/s]'); xlabel('t [s]');
legend('F_{1in}', 'F_D', 'Location', 'best');
exportgraphics(f, fullfile(folder, 'zad1_model_euler.png'), 'Resolution', 300); close(f);
