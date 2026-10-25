%%%% Linearyzacja modelu nieliniowego obiektu

clear; clc; close all;

%% Parametry
[model, F1_0, FD_0, h1_0, h2_0] = zad1_obiekt.parametry();
folder = fullfile(fileparts(mfilename('fullpath')), 'zad1_plots');

%% Linearyzacja w punkcie pracy
model_zlinearyzowany = zad1_obiekt.linearyzacja(model, h2_0);

%% Dyskretyzacja
Ts = 0.5;                 % krok [s]
Tk = 3000;                % czas symulacji [s]
N = round(Tk/Ts);
d = round(model.tau/Ts);     % opóźnienie w próbkach
t = (0:N)*Ts;

%% Eksperyment 1

F1in = zad1_scenariusze.przebieg_F1(t, F1_0);
FD = zad1_scenariusze.przebieg_FD(t, FD_0);
h2n = zad1_obiekt.sym_nieliniowy(model, F1in, FD, Ts, d, F1_0, h1_0, h2_0);
h2 = zad1_obiekt.sym_liniowy(model_zlinearyzowany, F1in, FD, Ts, d, F1_0, h1_0, h2_0);

f = zad1_wykresy.nowy();
subplot(2,1,1); plot(t, h2n, t, h2, '-'); grid on; ylabel('h_2 [cm]');
legend('Model nieliniowy', 'Model liniowy', 'Location', 'best');
subplot(2,1,2); stairs(t, F1in); hold on; stairs(t, FD, '-'); hold off;
grid on; ylabel('F [cm^3/s]'); xlabel('t [s]');
legend('F_{1in}', 'F_D', 'Location', 'best');
zad1_wykresy.zapisz(f, fullfile(folder, 'zad1_liniowy_vs_nieliniowy.png'));

%% Eksperyment 2

F1_dw = zad1_scenariusze.wzrost_spadek(t, F1_0);
FD_dw = FD_0*ones(size(t));
h2_dw = cell(1, 2);
for ik = 1:2
    h2_dw{ik} = {zad1_obiekt.sym_nieliniowy(model, F1_dw{ik}, FD_dw, Ts, d, F1_0, h1_0, h2_0), zad1_obiekt.sym_liniowy(model_zlinearyzowany, F1_dw{ik}, FD_dw, Ts, d, F1_0, h1_0, h2_0)};
end
zad1_wykresy.wzrost_spadek(t, F1_dw, F1_0, h2_dw, {'Model nieliniowy', 'Model liniowy'}, {'-', '-'}, fullfile(folder, 'zad1_liniowy_vs_nieliniowy_skoki.png'));

%% Eksperyment 3

skoki_F1 = [-0.5 -0.25 -0.1 0.1 0.25 0.5];   % względem F1_0
skoki_FD = [-0.5 0.5];                       % względem FD_0
[F1_s, FD_s, opis] = zad1_scenariusze.skoki_pojedyncze(t, F1_0, FD_0, skoki_F1, skoki_FD, 100);

nS  = numel(F1_s);
E   = zeros(1, nS);            % suma kwadratów różnicy h2 liniowy - nieliniowy
dh2 = zeros(2, nS);            % zmiana h2 w stanie ustalonym: [nieliniowy; liniowy]
for is = 1:nS
    h2n = zad1_obiekt.sym_nieliniowy(model, F1_s{is}, FD_s{is}, Ts, d, F1_0, h1_0, h2_0);
    h2  = zad1_obiekt.sym_liniowy(model_zlinearyzowany, F1_s{is}, FD_s{is}, Ts, d, F1_0, h1_0, h2_0);
    E(is) = sum((h2 - h2n).^2);
    dh2(:, is) = [h2n(end) - h2_0; h2(end) - h2_0];
end

%% Tabela jakości przybliżenia liniowego w zależności od wielkości skoku
fprintf('\nPorównanie modelu liniowego z nieliniowym (wyjście h2)\n');
fprintf('%10s %14s %14s %12s\n', 'skok', 'dh2 nieliniowy', 'dh2 liniowy', 'E');
for is = 1:nS
    fprintf('%10s %14.3f %14.3f %12.2f\n', opis{is}, dh2(1,is), dh2(2,is), E(is));
end
