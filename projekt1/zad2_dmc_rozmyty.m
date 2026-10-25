%%%% Rozmyty regulator DMC w wersji analitycznej z ograniczeniami.
%%%% Porównanie z konwencjonalnym DMC z zadania 1 oraz strojenie funkcji przynależności regulatora. 
%%%% Wymaga pliku zad2_model_rozmyty.mat tworzonego przez zad2_modele_rozmyte.m

clear; clc; close all;

%% Parametry obiektu
[model, F1_0, FD_0, h1_0, h2_0] = zad1_obiekt.parametry();

katalog = fileparts(mfilename('fullpath'));
folder = fullfile(katalog, 'zad2_plots');
plik_modelu = fullfile(katalog, 'zad2_model_rozmyty.mat');
if ~isfile(plik_modelu)
    error('Brak pliku %s - najpierw uruchom zad2_modele_rozmyte.m', plik_modelu);
end
dane = load(plik_modelu);

%% Parametry regulatorów
Tp_reg = 10;
limit_u = [0 110];
limit_du = [-10 10];

% Konwencjonalny DMC - nastawy jak w zadaniu 1
konw.D = 60; konw.N = 30; konw.Nu = 5; konw.lambda = 0.7;

% Rozmyty DMC - te same N, Nu, lambda; dłuższe D, bo przy wysokim h2 obiekt jest wolniejszy
rozm.D = 100; rozm.N = 30; rozm.Nu = 5; rozm.lambda = 0.7;

% Początkowe funkcje przynależności regulatora takie same jak w modelu obiektu
mf_start = dane.mf_popr;

% Warianty kształtu funkcji przynależności sprawdzane przy strojeniu regulatora
g_lista = [1 1.25 1.5 2 2.5];
szer_lista = [0.3 0.5 0.7 1];

%% Dyskretyzacja
Ts = 0.5;
Tk = 3000;
N_sim = round(Tk/Ts);
d = round(model.tau/Ts);
t = (0:N_sim)*Ts;
co = round(Tp_reg/Ts);
ip = 1:co:N_sim+1;

%% Scenariusz testowy
yzad = zad1_scenariusze.przebieg_yzad(t, h2_0);
FDc = zad1_scenariusze.przebieg_FD(t, FD_0);
scenariusz = zad1_dmc.utworz_scenariusz(model, yzad, FDc, Ts, co, d, F1_0, h1_0, h2_0, limit_u, limit_du);

%% Regulator konwencjonalny
regulator_konw = zad1_dmc.projekt(zad1_obiekt.linearyzacja(model, h2_0), konw, Ts, co, d);
[h2_k, u_k, E_k] = zad1_dmc.symulacja(regulator_konw, scenariusz);

%% Lokalne regulatory DMC zaprojektowane w centrach zbiorów
regulator_start = zad2_rozmyte.projekt_dmc(model, mf_start, rozm, Ts, co, d);
[h2_r, u_r, E_r] = zad1_dmc.symulacja(regulator_start, scenariusz);

fprintf('DMC konwencjonalny:                  E = %.2f\n', E_k);
fprintf('DMC rozmyty (%s): E = %.2f\n', zad2_rozmyte.opis(mf_start), E_r);

fprintf('\nOdpowiedzi skokowe modeli lokalnych po D = %d krokach:\n', rozm.D);
for i = 1:numel(mf_start.c)
    model_zlinearyzowany = zad1_obiekt.linearyzacja(model, mf_start.c(i));
    K = zad1_obiekt.wzmocnienie(model_zlinearyzowany);
    fprintf('h2 = %5.2f cm: K = %.4f, s(D) = %.2f%% K\n', model_zlinearyzowany.h2, K, 100*regulator_start.s(end,i)/K);
end


f = zad1_wykresy.nowy();
stairs((1:rozm.D)*Tp_reg, regulator_start.s); grid on; ylabel('s'); xlabel('t [s]');
legend(arrayfun(@(c) sprintf('model lokalny h_2 = %.2f cm', c), mf_start.c, 'UniformOutput', false), 'Location', 'best');
zad1_wykresy.zapisz(f, fullfile(folder, 'zad2_dmc_odp_skokowe_lokalne.png'));

nazwy = {sprintf('DMC konwencjonalny (E = %.1f)', E_k), sprintf('DMC rozmyty (E = %.1f)', E_r)};
zad1_wykresy.regulacja(t, ip, yzad, FDc, {h2_k, h2_r}, {u_k, u_r}, nazwy, nazwy, fullfile(folder, 'zad2_dmc_konw_vs_rozmyty.png'));

%% Strojenie funkcji przynależności regulatora
warianty = zad2_rozmyte.warianty(numel(mf_start.c), g_lista, szer_lista, dane.h2_min, dane.h2_max);
nW = numel(warianty);
E_war = zeros(1, nW);
for iw = 1:nW
    regulator = zad2_rozmyte.projekt_dmc(model, warianty{iw}, rozm, Ts, co, d);
    [~, ~, E_war(iw)] = zad1_dmc.symulacja(regulator, scenariusz);
end
[E_best, iw_best] = min(E_war);
if E_best < E_r
    mf_stroj = warianty{iw_best};
else
    mf_stroj = mf_start;
end
regulator_stroj = zad2_rozmyte.projekt_dmc(model, mf_stroj, rozm, Ts, co, d);
[h2_s, u_s, E_s] = zad1_dmc.symulacja(regulator_stroj, scenariusz);

fprintf('\nStrojenie funkcji przynależności regulatora rozmytego\n');
for iw = 1:nW
    fprintf('%-40s E = %10.2f\n', zad2_rozmyte.opis(warianty{iw}), E_war(iw));
end
fprintf('Przed strojeniem: %-40s E = %.2f\n', zad2_rozmyte.opis(mf_start), E_r);
fprintf('Po strojeniu:     %-40s E = %.2f\n', zad2_rozmyte.opis(mf_stroj), E_s);
fprintf('DMC konwencjonalny:                                        E = %.2f\n', E_k);

% Funkcje przynależności regulatora przed i po strojeniu
hh = linspace(dane.h2_min - 1, dane.h2_max + 1, 500);
zad2_rozmyte.wykres_mf(hh, {mf_start, mf_stroj}, {['Przed: ' zad2_rozmyte.opis(mf_start)], ['Po: ' zad2_rozmyte.opis(mf_stroj)]}, [1 2], fullfile(folder, 'zad2_dmc_funkcje_przynaleznosci.png'));

nazwy = {sprintf('DMC konwencjonalny (E = %.1f)', E_k), sprintf('DMC rozmyty przed strojeniem (E = %.1f)', E_r), sprintf('DMC rozmyty po strojeniu (E = %.1f)', E_s)};
zad1_wykresy.regulacja(t, ip, yzad, FDc, {h2_k, h2_r, h2_s}, {u_k, u_r, u_s}, nazwy, nazwy, fullfile(folder, 'zad2_dmc_przed_po_strojeniu.png'));
