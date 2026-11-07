%%%% Regulacja DMC zaprojektowana na modelu liniowym sterująca obiektem nieliniowym

clear; clc; close all;

%% Parametry obiektu
[model, F1_0, FD_0, h1_0, h2_0] = zad1_obiekt.parametry();

% Folder na wykresy
folder = fullfile(fileparts(mfilename('fullpath')), 'zad1_plots');

Tp_reg   = 10;
limit_u  = [0 110];
limit_du = [-10 10];

%% Nastawy bazowe i parametry do badania wpływu nastaw (D dobierane z odpowiedzi skokowej)
nastawy_bazowe.N      = 30;
nastawy_bazowe.Nu     = 5;
nastawy_bazowe.lambda = 0.7;
lambda_lista = [0.1 0.7 2 5];
N_lista = [10 20 30 60];
Nu_lista = [1 2 5 10];

%% Dyskretyzacja
Ts    = 0.5;
Tk    = 3000;
N_sim = round(Tk/Ts);
d     = round(model.tau/Ts);
t     = (0:N_sim)*Ts;
co    = round(Tp_reg/Ts);
ip    = 1:co:N_sim+1;

%% Model liniowy w punkcie pracy oraz scenariusz testowy
model_zlinearyzowany = zad1_obiekt.linearyzacja(model, h2_0);

yzad = zad1_scenariusze.przebieg_yzad(t, h2_0);
FDc  = zad1_scenariusze.przebieg_FD(t, FD_0);
scenariusz = zad1_dmc.utworz_scenariusz(model, yzad, FDc, Ts, co, d, F1_0, h1_0, h2_0, limit_u, limit_du);

%% Odpowiedź skokowa modelu liniowego - podstawa doboru horyzontu dynamiki D
D_max = 90;   % liczba próbek odpowiedzi skokowej na wykresie
s = zad1_dmc.odp_skokowa(model_zlinearyzowany, D_max, Ts, co, d);

% Wzmocnienie statyczne modelu liniowego
K_stat = zad1_obiekt.wzmocnienie(model_zlinearyzowany);

fprintf('Wzmocnienie statyczne zlinearyzowanego modelu: %.4f\n', K_stat);

f = zad1_wykresy.nowy();
stairs(1:D_max, s); hold on;
yline(K_stat, '--');
hold off;
grid on; ylabel('s'); xlabel('k');
title('Odpowiedź skokowa modelu liniowego', 'FontWeight', 'normal')
legend('Odpowiedź skokowa s', 'Wartość ustalona', 'Location', 'best');
zad1_wykresy.zapisz(f, fullfile(folder, 'zad1_dmc_odp_skokowa.png'));

% Horyzont dynamiki D dobrany na podstawie wykresu odpowiedzi skokowej
D = 80;
nastawy_bazowe.D = D;

%% Wpływ nastaw na jakość regulacji - zmiana jednego parametru, reszta bez zmian
pola     = {'N', 'Nu', 'lambda'};
etykiety = {'N', 'N_u', '\lambda'};
listy    = {N_lista, Nu_lista, lambda_lista};

for ib = 1:numel(pola)
    wart = listy{ib};
    h2b  = cell(1, numel(wart));
    ub   = cell(1, numel(wart));
    E_b  = zeros(size(wart));
    for iw = 1:numel(wart)
        nastawy_b = nastawy_bazowe;
        nastawy_b.(pola{ib}) = wart(iw);
        regulator_b = zad1_dmc.projekt(model_zlinearyzowany, nastawy_b, Ts, co, d);
        [h2b{iw}, ub{iw}, E_b(iw)] = zad1_dmc.symulacja(regulator_b, scenariusz);
    end
    etyk = arrayfun(@(w) sprintf('%s = %g', etykiety{ib}, w), wart, 'UniformOutput', false);
    zad1_wykresy.regulacja(t, ip, yzad, [], h2b, ub, etyk, {}, fullfile(folder, sprintf('zad1_dmc_strojenie_%s.png', pola{ib})), sprintf('Wpływ parametru %s na jakość regulacji DMC', etykiety{ib}));
    fprintf('\nWpływ %s na wskaźnik jakości E:\n', pola{ib});
    fprintf('%10g  E = %.2f\n', [wart; E_b]);
end

%% Wybrane najlepsze parametry regulatora DMC i regulacja obiektu nieliniowego
nastawy.D      = D;
nastawy.N      = 30;
nastawy.Nu     = 5;
nastawy.lambda = 0.1;

regulator = zad1_dmc.projekt(model_zlinearyzowany, nastawy, Ts, co, d);
[h2c, u, E] = zad1_dmc.symulacja(regulator, scenariusz);
fprintf('Wskaźnik jakości regulacji DMC: %.2f\n', E);

zad1_wykresy.regulacja(t, ip, yzad, FDc, {h2c}, {u}, {'Wyjście h_2'}, {'Sterowanie F_{1in}'}, fullfile(folder, 'zad1_dmc_regulacja.png'), 'Regulacja DMC obiektu nieliniowego');
