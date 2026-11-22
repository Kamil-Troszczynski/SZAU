%%%% Zadanie 4: numeryczny rozmyty regulator DMC typu NPL z ograniczeniami sterowania
%%%% i porównanie z regulatorami z poprzednich zadań:
%%%%   - DMC konwencjonalny (zad. 1), - rozmyty DMC analityczny (zad. 2), - rozmyty DMC-SL (zad. 3).
%%%% Wymaga pliku zad2_model_rozmyty.mat (zad2_modele_rozmyte.m) oraz quadprog (Optimization Toolbox).
%%%%
%%%% Trzy przypadki:
%%%%   1. skoki wartości zadanej, nastawy z zad. 2/3, ograniczenie przyrostu nieaktywne,
%%%%   2. szeroki zakres: skoki wartości zadanej +/-50% oraz skoki zakłócenia FD (jak w zad. 1/2),
%%%%   3. regulatory agresywne, ciasne ograniczenie przyrostu i górne ograniczenie wartości sterowania.

clear; clc; close all;

%% Parametry obiektu
[model, F1_0, FD_0, h1_0, h2_0] = zad1_obiekt.parametry();

katalog = fileparts(mfilename('fullpath'));
folder = fullfile(katalog, 'zad4_plots');
if ~isfolder(folder), mkdir(folder); end
plik_modelu = fullfile(katalog, 'zad2_model_rozmyty.mat');
if ~isfile(plik_modelu)
    error('Brak pliku %s - najpierw uruchom zad2_modele_rozmyte.m', plik_modelu);
end
dane = load(plik_modelu);
mf = dane.mf_popr;              % funkcje przynależności jak w modelu obiektu (PDC)

%% Parametry regulatorów
Tp_reg = 10;

% Konwencjonalny DMC - nastawy z zadania 1
konw.D = 80; konw.N = 30; konw.Nu = 5; konw.lambda = 0.1;

% Regulatory rozmyte - nastawy z zadań 2 i 3
rozm.D = 100; rozm.N = 30; rozm.Nu = 5; rozm.lambda = 0.1;

% Przypadek 3: bardzo mała kara za przyrosty, ciasne ograniczenia
lambda_3 = 0.0003;
du_max_3 = 2;
u_max_3  = 75;

typy  = {'konwencjonalny', 'analityczny', 'sl', 'npl'};
nazwy = {'DMC konwencjonalny', 'FDMC analityczny', 'FDMC-SL', 'FDMC-NPL'};

%% Dyskretyzacja
Ts = 0.5;
d = round(model.tau/Ts);
co = round(Tp_reg/Ts);

%% Przypadek 1: skoki wartości zadanej, stałe zakłócenie
t1 = 0:Ts:1550;
yzad1 = zad1_scenariusze.przebieg(t1, h2_0, [50 550 1050], [1.3 1.0 0.7]);
sc1 = zad1_dmc.utworz_scenariusz(model, yzad1, FD_0*ones(size(t1)), Ts, co, d, F1_0, h1_0, h2_0, [0 110], [-10 10]);

fprintf('=== Przypadek 1: skoki wartości zadanej, lambda = %g, |du| <= %g ===\n', rozm.lambda, 10);
wyn1 = symuluj(typy, nazwy, model, mf, konw, rozm, Ts, co, d, F1_0, FD_0, h2_0, sc1);
drukuj(wyn1);
wykres_przypadku(t1, co, sc1, wyn1, fullfile(folder, 'zad4_przypadek1_skoki_yzad.png'), ...
    sprintf('Skoki wartości zadanej (\\lambda = %g, |\\Deltau| \\leq %g)', rozm.lambda, 10));

%% Przypadek 2: szeroki zakres wartości zadanej (+/-50%) i skoki zakłócenia FD
t2 = 0:Ts:3000;
yzad2 = zad1_scenariusze.przebieg_yzad(t2, h2_0);
FD2 = zad1_scenariusze.przebieg_FD(t2, FD_0);
sc2 = zad1_dmc.utworz_scenariusz(model, yzad2, FD2, Ts, co, d, F1_0, h1_0, h2_0, [0 110], [-10 10]);

fprintf('\n=== Przypadek 2: wartość zadana +/-50%%, skoki FD, lambda = %g, |du| <= %g ===\n', rozm.lambda, 10);
wyn2 = symuluj(typy, nazwy, model, mf, konw, rozm, Ts, co, d, F1_0, FD_0, h2_0, sc2);
drukuj(wyn2);
wykres_przypadku(t2, co, sc2, wyn2, fullfile(folder, 'zad4_przypadek2_szeroki_zakres_zaklocenia.png'), ...
    sprintf('Wartość zadana \\pm50%% i zakłócenie F_D (\\lambda = %g, |\\Deltau| \\leq %g)', rozm.lambda, 10));

%% Przypadek 3: regulatory agresywne, ograniczenia przyrostu i wartości sterowania aktywne
sc3 = sc1;
sc3.limit_u  = [0 u_max_3];
sc3.limit_du = du_max_3*[-1 1];
konw3 = konw; konw3.lambda = lambda_3;
rozm3 = rozm; rozm3.lambda = lambda_3;

fprintf('\n=== Przypadek 3: lambda = %g, |du| <= %g, u <= %g ===\n', lambda_3, du_max_3, u_max_3);
wyn3 = symuluj(typy, nazwy, model, mf, konw3, rozm3, Ts, co, d, F1_0, FD_0, h2_0, sc3);
drukuj(wyn3);
wykres_przypadku(t1, co, sc3, wyn3, fullfile(folder, 'zad4_przypadek3_ograniczenia_aktywne.png'), ...
    sprintf('Ograniczenia aktywne (\\lambda = %g, |\\Deltau| \\leq %g, u \\leq %g)', lambda_3, du_max_3, u_max_3));

%% Funkcje pomocnicze
% Symulacja wszystkich regulatorów w jednym scenariuszu; wynik - tablica struktur
function wyn = symuluj(typy, nazwy, model, mf, konw, rozm, Ts, co, d, F1_0, FD_0, h2_0, scenariusz)
    for i = numel(typy):-1:1
        if strcmp(typy{i}, 'konwencjonalny')
            regulator = zad1_dmc.projekt(zad1_obiekt.linearyzacja(model, h2_0), konw, Ts, co, d);
            regulator.typ = 'analityczny';
        else
            regulator = zad4_npl.projekt(model, mf, rozm, Ts, co, d, typy{i}, F1_0, FD_0);
        end
        [h2, u, E, czas] = zad4_npl.symulacja(regulator, scenariusz);
        [TV, n_ogr] = zad3_sl.wskazniki(u, scenariusz);
        przereg = zad3_sl.przeregulowanie(h2, scenariusz.yzad);
        wyn(i) = struct('nazwa', nazwy{i}, 'h2', h2, 'u', u, 'E', E, 'przereg', przereg, ...
                        'TV', TV, 'n_ogr', n_ogr, 'czas', czas);
    end
end

function drukuj(wyn)
    for w = wyn
        fprintf('  %-20s E = %8.2f   przereg. = %4.2f cm   TV = %6.1f   kroki na ogr. = %3d   czas = %7.3f ms\n', ...
            w.nazwa, w.E, w.przereg, w.TV, w.n_ogr, w.czas);
    end
end

% Przebiegi regulatorów: wyjście, sterowanie, przyrost sterowania z ograniczeniem
% oraz (gdy się zmienia) zakłócenie FD
function wykres_przypadku(t, co, sc, wyn, plik, tytul)
    ip = 1:co:numel(t);
    tt = t(ip);
    kolory = [0.5 0.5 0.5; 0 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188];
    style = {{':', 'LineWidth', 1.5}, {'-', 'LineWidth', 1.5}, {'--', 'LineWidth', 1.5}, {'-', 'LineWidth', 2}};
    rysuj = @(x, y, i) stairs(x, y, style{i}{:}, 'Color', kolory(i,:));
    z_FD = any(sc.FD ~= sc.FD(1));

    f = zad1_wykresy.nowy(1000, 650 + 150*z_FD);
    uklad = tiledlayout(f, 3 + z_FD, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(uklad, tytul, 'FontSize', 12);

    ax(1) = nexttile; hold on;
    for i = 1:numel(wyn), rysuj(tt, wyn(i).h2(ip), i); end
    stairs(tt, sc.yzad(ip), 'k-.', 'LineWidth', 1);
    ylabel('h_2 [cm]');
    lgd = legend([{wyn.nazwa}, {'Wartość zadana'}], 'Orientation', 'horizontal');
    lgd.Layout.Tile = 'north';

    ax(2) = nexttile; hold on;
    for i = 1:numel(wyn), rysuj(tt, wyn(i).u(ip), i); end
    if sc.limit_u(2) < 110
        yline(sc.limit_u(2), 'k-.', 'LineWidth', 1);
    end
    ylabel('F_{1in} [cm^3/s]');

    ax(3) = nexttile; hold on;
    du_max = sc.limit_du(2);
    for i = 1:numel(wyn), rysuj(tt(2:end), diff(wyn(i).u(ip)), i); end
    yline([-du_max du_max], 'k-.', 'LineWidth', 1);
    ylim(1.35*du_max*[-1 1]);
    ylabel('\Deltau [cm^3/s]');

    if z_FD
        ax(4) = nexttile;
        stairs(tt, sc.FD(ip), 'k', 'LineWidth', 1.2);
        ylabel('F_D [cm^3/s]');
    end
    xlabel(ax(end), 't [s]');

    for a = ax
        grid(a, 'on'); box(a, 'on'); xlim(a, [tt(1) tt(end)]);
    end
    set(ax(1:end-1), 'XTickLabel', []);
    linkaxes(ax, 'x');
    zad1_wykresy.zapisz(f, plik);
end
