%%%% Model obiektu rozpisany metodą Eulera (krok metody: zad1_obiekt.krok)

clear; clc; close all;

%% Parametry
[model, F1_0, FD_0, h1_0, h2_0] = zad1_obiekt.parametry();
folder = fullfile(fileparts(mfilename('fullpath')), 'zad1_plots');

%% Dyskretyzacja
Ts = 0.5;                 % krok [s]
Tk = 3000;                % czas symulacji [s]
N = round(Tk/Ts);
d = round(model.tau/Ts);     % opóźnienie w próbkach
t = (0:N)*Ts;

%% Wejścia
F1in = zad1_scenariusze.przebieg_F1(t, F1_0);
FD = zad1_scenariusze.przebieg_FD(t, FD_0);

%% Symulacja
h2 = zad1_obiekt.sym_nieliniowy(model, F1in, FD, Ts, d, F1_0, h1_0, h2_0);

%% Wykresy
f = zad1_wykresy.nowy();
subplot(2,1,1); plot(t, h2); grid on; ylabel('h_2 [cm]');
subplot(2,1,2); stairs(t, F1in); hold on; stairs(t, FD, 'black'); hold off;
grid on; ylabel('F [cm^3/s]'); xlabel('t [s]');
legend('F_{1in}', 'F_D', 'Location', 'best');
zad1_wykresy.zapisz(f, fullfile(folder, 'zad1_model_euler.png'));
