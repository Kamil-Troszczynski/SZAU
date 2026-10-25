%%%% Modele rozmyte Takagi-Sugeno (2, 3, 4, 5 modeli lokalnych), porównanie z modelem
%%%% liniowym i nieliniowym oraz poprawa wybranego modelu przez zmianę kształtu funkcji
%%%% przynależności.

clear; clc; close all;

%% Parametry obiektu
[model, F1_0, FD_0, h1_0, h2_0] = zad1_obiekt.parametry();

%% Parametry modeli rozmytych
L_lista = [2 3 4 5];        % liczba modeli lokalnych
L_wyb = 5;                  % model wybrany do dalszych badań (wybrany na podstawie tabeli "Błąd modeli względem modelu nieliniowego (suma kwadratów różnic h2"))

% Zmienna harmonogramująca
h2_min = ((0.5*F1_0 + FD_0)/model.a2)^2;
h2_max = ((1.5*F1_0 + FD_0)/model.a2)^2;

% Warianty kształtu funkcji przynależności sprawdzane przy poprawie modelu
g_lista = [1 1.25 1.5 2 2.5];    % rozkład centrów: 1 - równomierny, > 1 - gęściej przy małych h2
szer_lista = [0.3 0.5 0.7 1];       % szerokość funkcji Gaussa względem odstępu centrów

katalog = fileparts(mfilename('fullpath'));
folder = fullfile(katalog, 'zad2_plots');
plik_modelu = fullfile(katalog, 'zad2_model_rozmyty.mat');

if ~any(L_lista == L_wyb)
    error('L_wyb = %d musi należeć do L_lista', L_wyb);
end

%% Dyskretyzacja
Ts = 0.5;
Tk = 3000;
N = round(Tk/Ts);
d = round(model.tau/Ts);
t = (0:N)*Ts;

%% Model liniowy to model Takagi - Sugeko z jednym modelem lokalnym w punkcie pracy
model_zlinearyzowany = zad1_obiekt.linearyzacja(model, h2_0);
mf_lin = zad2_rozmyte.mf_liniowy(h2_0);

%% Modele rozmyte: równomierny podział h2, trójkątne funkcje przynależności
nL = numel(L_lista);
mf_L = arrayfun(@(L) zad2_rozmyte.utworz_mf(L, 'trojkatne', 1, 0, h2_min, h2_max), L_lista, 'UniformOutput', false);

hh = linspace(h2_min - 1, h2_max + 1, 500);
zad2_rozmyte.wykres_mf(hh, mf_L, arrayfun(@(L) sprintf('%d modele lokalne', L), L_lista, 'UniformOutput', false), [2 2], fullfile(folder, 'zad2_funkcje_przynaleznosci.png'));

%% Scenariusze testowe
% a) pojedyncze skoki z punktu pracy - do tabeli błędów i do poprawy modelu
skoki_F1 = [-0.5 -0.25 -0.1 0.1 0.25 0.5];   % względem F1_0
skoki_FD = [-0.5 0.5];                       % względem FD_0
[F1_s, FD_s, opis_s] = zad1_scenariusze.skoki_pojedyncze(t, F1_0, FD_0, skoki_F1, skoki_FD, 100);
nS = numel(F1_s);

% b) skoki F1in w dwóch miejscach, wzrost i spadek - do wykresów (jak w zadaniu 1)
F1_dw = zad1_scenariusze.wzrost_spadek(t, F1_0);
FD_dw = FD_0*ones(size(t));

% Odpowiedzi modelu nieliniowego - punkt odniesienia
h2n_s = cellfun(@(F1, FD) zad1_obiekt.sym_nieliniowy(model, F1, FD, Ts, d, F1_0, h1_0, h2_0), F1_s, FD_s, 'UniformOutput', false);
h2n_dw = cellfun(@(F1) zad1_obiekt.sym_nieliniowy(model, F1, FD_dw, Ts, d, F1_0, h1_0, h2_0), F1_dw, 'UniformOutput', false);

%% Porównanie modeli rozmytych z modelem liniowym i nieliniowym
nM = nL + 1;
modele_lok_lista = [{model_zlinearyzowany}, cellfun(@(mf) zad2_rozmyte.modele_lokalne(model, mf.c), mf_L, 'UniformOutput', false)];
mf_all = [{mf_lin}, mf_L];
nazwy = [{'liniowy'}, arrayfun(@(L) sprintf('TS%d', L), L_lista, 'UniformOutput', false)];

E = zeros(nS, nM);
for im = 1:nM
    for is = 1:nS
        h2m = zad2_rozmyte.sym_ts(modele_lok_lista{im}, mf_all{im}, F1_s{is}, FD_s{is}, Ts, d, F1_0, h1_0, h2_0);
        E(is, im) = sum((h2m - h2n_s{is}).^2);
    end
end

fprintf('\nBłąd modeli względem modelu nieliniowego (suma kwadratów różnic h2)\n');
fprintf('%-10s', 'skok'); fprintf('%12s', nazwy{:}); fprintf('\n');
for is = 1:nS
    fprintf('%-10s', opis_s{is}); fprintf('%12.2f', E(is,:)); fprintf('\n');
end
fprintf('%-10s', 'suma'); fprintf('%12.2f', sum(E, 1)); fprintf('\n');

h2_dw = cell(1, 2);
for ik = 1:2
    h2_dw{ik} = [h2n_dw(ik), cellfun(@(modele_lok, mf) zad2_rozmyte.sym_ts(modele_lok, mf, F1_dw{ik}, FD_dw, Ts, d, F1_0, h1_0, h2_0), modele_lok_lista, mf_all, 'UniformOutput', false)];
end
zad1_wykresy.wzrost_spadek(t, F1_dw, F1_0, h2_dw, [{'nieliniowy'}, nazwy], [{'k', '--'}, repmat({''}, 1, nL)], fullfile(folder, 'zad2_porownanie_modeli.png'));

%% Poprawa wybranego modelu przez zmianę kształtu funkcji przynależności
mf_wyb = mf_L{L_lista == L_wyb};
warianty = zad2_rozmyte.warianty(L_wyb, g_lista, szer_lista, h2_min, h2_max);
nW = numel(warianty);

E_war = zeros(1, nW);
for iw = 1:nW
    modele_lok = zad2_rozmyte.modele_lokalne(model, warianty{iw}.c);
    for is = 1:nS
        h2m = zad2_rozmyte.sym_ts(modele_lok, warianty{iw}, F1_s{is}, FD_s{is}, Ts, d, F1_0, h1_0, h2_0);
        E_war(iw) = E_war(iw) + sum((h2m - h2n_s{is}).^2);
    end
end
[~, iw_best] = min(E_war);
mf_popr = warianty{iw_best};

E_przed = sum(E(:, find(L_lista == L_wyb) + 1));
fprintf('\nPoprawa modelu TS%d - suma błędów po wszystkich skokach\n', L_wyb);
for iw = 1:nW
    fprintf('%-40s E = %10.2f\n', zad2_rozmyte.opis(warianty{iw}), E_war(iw));
end
fprintf('Przed poprawą: %-40s E = %.2f\n', zad2_rozmyte.opis(mf_wyb), E_przed);
fprintf('Po poprawie:   %-40s E = %.2f\n', zad2_rozmyte.opis(mf_popr), E_war(iw_best));

% Funkcje przynależności przed i po poprawie
zad2_rozmyte.wykres_mf(hh, {mf_wyb, mf_popr}, {['Przed: ' zad2_rozmyte.opis(mf_wyb)], ['Po: ' zad2_rozmyte.opis(mf_popr)]}, [1 2], fullfile(folder, 'zad2_funkcje_przynaleznosci_poprawa.png'));

% Odpowiedzi modelu przed i po poprawie
modele_lok_wyb  = zad2_rozmyte.modele_lokalne(model, mf_wyb.c);
modele_lok_popr = zad2_rozmyte.modele_lokalne(model, mf_popr.c);
for ik = 1:2
    h2_dw{ik} = {h2n_dw{ik}, ...
        zad2_rozmyte.sym_ts(model_zlinearyzowany,  mf_lin,  F1_dw{ik}, FD_dw, Ts, d, F1_0, h1_0, h2_0), ...
        zad2_rozmyte.sym_ts(modele_lok_wyb,  mf_wyb,  F1_dw{ik}, FD_dw, Ts, d, F1_0, h1_0, h2_0), ...
        zad2_rozmyte.sym_ts(modele_lok_popr, mf_popr, F1_dw{ik}, FD_dw, Ts, d, F1_0, h1_0, h2_0)};
end
zad1_wykresy.wzrost_spadek(t, F1_dw, F1_0, h2_dw, {'nieliniowy', 'liniowy', sprintf('TS%d przed poprawą', L_wyb), sprintf('TS%d po poprawie', L_wyb)}, {'k', '--', '', ''}, fullfile(folder, 'zad2_model_przed_po_poprawie.png'));

%% Zapis funkcji przynależności wybranego modelu
save(plik_modelu, 'mf_wyb', 'mf_popr', 'h2_min', 'h2_max');
fprintf('\nZapisano funkcje przynależności do %s\n', plik_modelu);
