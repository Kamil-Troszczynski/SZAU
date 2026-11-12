%%%% Zadanie 3: numeryczny rozmyty regulator DMC typu SL z ograniczeniami sterowania
%%%% i porównanie z rozmytym regulatorem DMC analitycznym z zadania 2.
%%%% Wymaga pliku zad2_model_rozmyty.mat (zad2_modele_rozmyte.m) oraz quadprog (Optimization Toolbox).
%%%%
%%%% Dwa przypadki dla tego samego przebiegu wartości zadanej:
%%%%   1. ograniczenie przyrostu sterowania nie jest osiągane (nastawy z zad. 2),
%%%%   2. regulator agresywny i ciasne ograniczenie przyrostu - ograniczenie silnie wpływa na sterowanie.

clear; clc; close all;

%% Parametry obiektu
[model, F1_0, FD_0, h1_0, h2_0] = zad1_obiekt.parametry();

katalog = fileparts(mfilename('fullpath'));
folder = fullfile(katalog, 'zad3_plots');
if ~isfolder(folder), mkdir(folder); end
plik_modelu = fullfile(katalog, 'zad2_model_rozmyty.mat');
if ~isfile(plik_modelu)
    error('Brak pliku %s - najpierw uruchom zad2_modele_rozmyte.m', plik_modelu);
end
dane = load(plik_modelu);
mf = dane.mf_popr;              % funkcje przynależności jak w modelu obiektu (PDC)

%% Parametry regulatorów
Tp_reg = 10;
limit_u = [0 110];

% Przypadek 1: nastawy jak dla rozmytego DMC w zadaniu 2
nastawy_1.D = 100; nastawy_1.N = 30; nastawy_1.Nu = 5; nastawy_1.lambda = 0.1;
du_max_1 = 10;

% Przypadek 2: bardzo mała kara za przyrosty sterowania i ciasne ograniczenie przyrostu
nastawy_2 = nastawy_1;
nastawy_2.lambda = 0.0003;
du_max_2 = 2;

typy  = {'analityczny', 'sl'};
nazwy = {'DMC analityczny (przycięcie)', 'DMC SL numeryczny (QP)'};

%% Dyskretyzacja i scenariusz: trzy skoki wartości zadanej, stałe zakłócenie
Ts = 0.5;
d = round(model.tau/Ts);
co = round(Tp_reg/Ts);
t = 0:Ts:1550;
yzad = zad1_scenariusze.przebieg(t, h2_0, [50 550 1050], [1.3 1.0 0.7]);
scenariusz = zad1_dmc.utworz_scenariusz(model, yzad, FD_0*ones(size(t)), Ts, co, d, F1_0, h1_0, h2_0, limit_u, du_max_1*[-1 1]);

%% Przypadek 1: ograniczenie nie wpływa na sterowanie
fprintf('=== Przypadek 1: lambda = %g, |du| <= %g ===\n', nastawy_1.lambda, du_max_1);
wyn = symuluj(typy, nazwy, model, mf, nastawy_1, Ts, co, d, scenariusz);
drukuj(wyn);
wykres_przypadku(t, co, yzad, wyn, du_max_1, fullfile(folder, 'zad3_przypadek1_ograniczenie_nieaktywne.png'), ...
    sprintf('Ograniczenie nie wpływa na sterowanie (\\lambda = %g, |\\Deltau| \\leq %g)', nastawy_1.lambda, du_max_1));

%% Przypadek 2: ograniczenie silnie wpływa na sterowanie
scenariusz.limit_du = du_max_2*[-1 1];
fprintf('\n=== Przypadek 2: lambda = %g, |du| <= %g ===\n', nastawy_2.lambda, du_max_2);
wyn = symuluj(typy, nazwy, model, mf, nastawy_2, Ts, co, d, scenariusz);
drukuj(wyn);
wykres_przypadku(t, co, yzad, wyn, du_max_2, fullfile(folder, 'zad3_przypadek2_ograniczenie_aktywne.png'), ...
    sprintf('Ograniczenie silnie wpływa na sterowanie (\\lambda = %g, |\\Deltau| \\leq %g)', nastawy_2.lambda, du_max_2));

%% Funkcje pomocnicze
% Symulacja kilku regulatorów w jednym scenariuszu; wynik - tablica struktur
function wyn = symuluj(typy, nazwy, model, mf, nastawy, Ts, co, d, scenariusz)
    for i = numel(typy):-1:1
        regulator = zad3_sl.projekt(model, mf, nastawy, Ts, co, d, typy{i});
        [h2, u, E, czas] = zad3_sl.symulacja(regulator, scenariusz);
        [TV, n_ogr] = zad3_sl.wskazniki(u, scenariusz);
        przereg = zad3_sl.przeregulowanie(h2, scenariusz.yzad);
        wyn(i) = struct('nazwa', nazwy{i}, 'h2', h2, 'u', u, 'E', E, 'przereg', przereg, ...
                        'TV', TV, 'n_ogr', n_ogr, 'czas', czas);
    end
end

function drukuj(wyn)
    for w = wyn
        fprintf('  %-30s E = %7.2f   przereg. = %4.2f cm   TV = %5.1f   kroki na ogr. = %2d   czas = %5.3f ms\n', ...
            w.nazwa, w.E, w.przereg, w.TV, w.n_ogr, w.czas);
    end
end

% Przebiegi dwóch regulatorów: wyjście, sterowanie i przyrost sterowania z ograniczeniem
function wykres_przypadku(t, co, yzad, wyn, du_max, plik, tytul)
    ip = 1:co:numel(t);
    tt = t(ip);
    kolory = [0 0.447 0.741; 0.850 0.325 0.098];
    style = {{'-', 'LineWidth', 2}, {'--', 'LineWidth', 2}};
    rysuj = @(x, y, i) stairs(x, y, style{i}{:}, 'Color', kolory(i,:));

    f = zad1_wykresy.nowy(1000, 650);
    uklad = tiledlayout(f, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(uklad, tytul, 'FontSize', 12);

    ax(1) = nexttile; hold on;
    for i = 1:numel(wyn), rysuj(tt, wyn(i).h2(ip), i); end
    stairs(tt, yzad(ip), 'k:', 'LineWidth', 1.2);
    ylabel('h_2 [cm]');
    lgd = legend([{wyn.nazwa}, {'Wartość zadana'}], 'Orientation', 'horizontal');
    lgd.Layout.Tile = 'north';

    ax(2) = nexttile; hold on;
    for i = 1:numel(wyn), rysuj(tt, wyn(i).u(ip), i); end
    ylabel('u [cm^3/s]');

    ax(3) = nexttile; hold on;
    for i = 1:numel(wyn), rysuj(tt(2:end), diff(wyn(i).u(ip)), i); end
    yline([-du_max du_max], 'k-.', 'LineWidth', 1);
    text(tt(end), du_max, sprintf('|\\Deltau|_{max} = %g ', du_max), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 9);
    ylim(1.35*du_max*[-1 1]);
    ylabel('\Deltau [cm^3/s]'); xlabel('t [s]');

    for a = ax
        grid(a, 'on'); box(a, 'on'); xlim(a, [tt(1) tt(end)]);
    end
    set(ax(1:2), 'XTickLabel', []);
    linkaxes(ax, 'x');
    zad1_wykresy.zapisz(f, plik);
end
