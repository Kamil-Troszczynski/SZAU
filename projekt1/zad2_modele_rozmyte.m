%%%% Modele rozmyte Takagi-Sugeno (2, 3, 4, 5 modeli lokalnych), porównanie z modelem
%%%% liniowym i nieliniowym oraz poprawa wybranego modelu przez zmianę kształtu funkcji
%%%% przynależności. Funkcje przynależności zapisywane do zad2_model_rozmyty.mat
%%%% (korzysta z nich zad2_dmc_rozmyty.m)
clear; clc; close all;

%% Parametry obiektu
p.A2 = 430; p.C1 = 0.85; p.a1 = 18; p.a2 = 24; p.tau = 70;
F1_0 = 55; FD_0 = 10;
h1_0 = ((F1_0 + FD_0)/p.a1)^2;
h2_0 = ((F1_0 + FD_0)/p.a2)^2;

%% Parametry modeli rozmytych
L_lista = [2 3 4 5];        % liczba modeli lokalnych
L_wyb   = 5;                % model wybrany do dalszych badań (wybór na podstawie tabeli błędów)

% Zmienna harmonogramująca: wyjście h2, zakres wynikający ze zmian F1 o +/-50% od punktu pracy
h2_min = ((0.5*F1_0 + FD_0)/p.a2)^2;
h2_max = ((1.5*F1_0 + FD_0)/p.a2)^2;

% Warianty kształtu funkcji przynależności sprawdzane przy poprawie modelu
g_lista    = [1 1.25 1.5 2 2.5];    % rozkład centrów: 1 - równomierny, > 1 - gęściej przy małych h2
szer_lista = [0.3 0.5 0.7 1];       % szerokość funkcji Gaussa względem odstępu centrów

% Folder na wykresy i plik z wynikiem
katalog     = fileparts(mfilename('fullpath'));
folder      = fullfile(katalog, 'zad2_plots');
plik_modelu = fullfile(katalog, 'zad2_model_rozmyty.mat');

if ~any(L_lista == L_wyb)
    error('L_wyb = %d musi należeć do L_lista', L_wyb);
end

%% Dyskretyzacja
Ts = 0.5;
Tk = 3000;
N  = round(Tk/Ts);
d  = round(p.tau/Ts);
t  = (0:N)*Ts;

%% Model liniowy = model TS z jednym modelem lokalnym w punkcie pracy
lok_lin = zad2_lokalny(p, h2_0);
mf_lin  = struct('L', 1, 'typ', 'trojkatne', 'g', 1, 'szer', 0, 'c', h2_0);

%% Modele rozmyte: równomierny podział h2, trójkątne funkcje przynależności
nL   = numel(L_lista);
mf_L = cell(1, nL);
for il = 1:nL
    mf_L{il} = zad2_mf(L_lista(il), 'trojkatne', 1, 0, h2_min, h2_max);
end

hh = linspace(h2_min - 1, h2_max + 1, 500);
f = figure('Color', 'w', 'Visible', 'off');
for il = 1:nL
    subplot(2,2,il); plot(hh, wartosci_mf(hh, mf_L{il})); grid on;
    xlabel('h_2 [cm]'); ylabel('\mu');
    title(sprintf('%d modele lokalne', L_lista(il)));
end
exportgraphics(f, fullfile(folder, 'zad2_funkcje_przynaleznosci.png'), 'Resolution', 300); close(f);

%% Scenariusze testowe
% a) pojedyncze skoki z punktu pracy - do tabeli błędów i do poprawy modelu
skoki_F1 = [-0.5 -0.25 -0.1 0.1 0.25 0.5];   % względem F1_0
skoki_FD = [-0.5 0.5];                       % względem FD_0
t_skok   = 100;

nS     = numel(skoki_F1) + numel(skoki_FD);
F1_s   = cell(1, nS);
FD_s   = cell(1, nS);
opis_s = cell(1, nS);
for is = 1:nS
    F1_s{is} = F1_0*ones(1, N+1);
    FD_s{is} = FD_0*ones(1, N+1);
    if is <= numel(skoki_F1)
        F1_s{is}(t >= t_skok) = (1 + skoki_F1(is))*F1_0;
        opis_s{is} = sprintf('F1 %+g%%', 100*skoki_F1(is));
    else
        sk = skoki_FD(is - numel(skoki_F1));
        FD_s{is}(t >= t_skok) = (1 + sk)*FD_0;
        opis_s{is} = sprintf('FD %+g%%', 100*sk);
    end
end

% b) skoki F1in w dwóch miejscach, wzrost i spadek - do wykresów (jak w zadaniu 1)
kierunki = [1 -1];
tytuly   = {'Wzrost F_{1in}', 'Spadek F_{1in}'};
FD_dw    = FD_0*ones(1, N+1);
F1_dw    = cell(1, 2);
for ik = 1:2
    F1_dw{ik} = F1_0*ones(1, N+1);
    F1_dw{ik}(t >= 100)  = (1 + kierunki(ik)*0.25)*F1_0;
    F1_dw{ik}(t >= 1500) = (1 + kierunki(ik)*0.5)*F1_0;
end

% Odpowiedzi modelu nieliniowego - punkt odniesienia
h2n_s = cell(1, nS);
for is = 1:nS
    h2n_s{is} = sym_nielin(p, F1_s{is}, FD_s{is}, Ts, d, F1_0, h1_0, h2_0);
end
h2n_dw = cell(1, 2);
for ik = 1:2
    h2n_dw{ik} = sym_nielin(p, F1_dw{ik}, FD_dw, Ts, d, F1_0, h1_0, h2_0);
end

%% Porównanie modeli rozmytych z modelem liniowym i nieliniowym
nM      = nL + 1;                  % model liniowy + modele TS
lok_all = cell(1, nM);
mf_all  = cell(1, nM);
lok_all{1} = lok_lin;
mf_all{1}  = mf_lin;
for il = 1:nL
    lok_all{il+1} = modele_lokalne(p, mf_L{il}.c);
    mf_all{il+1}  = mf_L{il};
end
nazwy = [{'liniowy'}, arrayfun(@(L) sprintf('TS%d', L), L_lista, 'UniformOutput', false)];

E = zeros(nS, nM);                 % suma kwadratów różnicy h2 względem modelu nieliniowego
for im = 1:nM
    for is = 1:nS
        h2m = sym_ts(lok_all{im}, mf_all{im}, F1_s{is}, FD_s{is}, Ts, d, F1_0, h1_0, h2_0);
        E(is, im) = sum((h2m - h2n_s{is}).^2);
    end
end

fprintf('\nBłąd modeli względem modelu nieliniowego (suma kwadratów różnic h2)\n');
fprintf('%-10s', 'skok'); fprintf('%12s', nazwy{:}); fprintf('\n');
for is = 1:nS
    fprintf('%-10s', opis_s{is}); fprintf('%12.2f', E(is,:)); fprintf('\n');
end
fprintf('%-10s', 'suma'); fprintf('%12.2f', sum(E, 1)); fprintf('\n');

f = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1200 700]);
for ik = 1:2
    subplot(2,2,ik); plot(t, h2n_dw{ik}, 'k'); hold on;
    for im = 1:nM
        h2m = sym_ts(lok_all{im}, mf_all{im}, F1_dw{ik}, FD_dw, Ts, d, F1_0, h1_0, h2_0);
        if im == 1
            plot(t, h2m, '--');
        else
            plot(t, h2m);
        end
    end
    hold off; grid on; ylabel('h_2 [cm]'); title(tytuly{ik});
    legend([{'nieliniowy'}, nazwy], 'Location', 'best');

    subplot(2,2,ik+2);
    stairs(t, F1_dw{ik}); grid on; ylabel('F_{1in} [cm^3/s]'); xlabel('t [s]');
    ylim([0.4 1.6]*F1_0);
end
exportgraphics(f, fullfile(folder, 'zad2_porownanie_modeli.png'), 'Resolution', 300); close(f);

%% Poprawa wybranego modelu przez zmianę kształtu funkcji przynależności
mf_wyb   = mf_L{L_lista == L_wyb};
warianty = zad2_warianty_mf(L_wyb, g_lista, szer_lista, h2_min, h2_max);
nW       = numel(warianty);

E_war = zeros(1, nW);              % suma błędów po wszystkich pojedynczych skokach
for iw = 1:nW
    lok = modele_lokalne(p, warianty{iw}.c);
    for is = 1:nS
        h2m = sym_ts(lok, warianty{iw}, F1_s{is}, FD_s{is}, Ts, d, F1_0, h1_0, h2_0);
        E_war(iw) = E_war(iw) + sum((h2m - h2n_s{is}).^2);
    end
end
[~, iw_best] = min(E_war);
mf_popr = warianty{iw_best};

E_przed = sum(E(:, find(L_lista == L_wyb) + 1));
fprintf('\nPoprawa modelu TS%d - suma błędów po wszystkich skokach\n', L_wyb);
for iw = 1:nW
    fprintf('%-40s E = %10.2f\n', zad2_opis_mf(warianty{iw}), E_war(iw));
end
fprintf('Przed poprawą: %-40s E = %.2f\n', zad2_opis_mf(mf_wyb), E_przed);
fprintf('Po poprawie:   %-40s E = %.2f\n', zad2_opis_mf(mf_popr), E_war(iw_best));

% Funkcje przynależności przed i po poprawie
f = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1100 400]);
subplot(1,2,1); plot(hh, wartosci_mf(hh, mf_wyb)); grid on;
xlabel('h_2 [cm]'); ylabel('\mu'); title(['Przed: ' zad2_opis_mf(mf_wyb)]);
subplot(1,2,2); plot(hh, wartosci_mf(hh, mf_popr)); grid on;
xlabel('h_2 [cm]'); ylabel('\mu'); title(['Po: ' zad2_opis_mf(mf_popr)]);
exportgraphics(f, fullfile(folder, 'zad2_funkcje_przynaleznosci_poprawa.png'), 'Resolution', 300); close(f);

% Odpowiedzi modelu przed i po poprawie
lok_wyb  = modele_lokalne(p, mf_wyb.c);
lok_popr = modele_lokalne(p, mf_popr.c);
f = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1200 700]);
for ik = 1:2
    h2_l  = sym_ts(lok_lin,  mf_lin,  F1_dw{ik}, FD_dw, Ts, d, F1_0, h1_0, h2_0);
    h2_pr = sym_ts(lok_wyb,  mf_wyb,  F1_dw{ik}, FD_dw, Ts, d, F1_0, h1_0, h2_0);
    h2_po = sym_ts(lok_popr, mf_popr, F1_dw{ik}, FD_dw, Ts, d, F1_0, h1_0, h2_0);

    subplot(2,2,ik);
    plot(t, h2n_dw{ik}, 'k', t, h2_l, '--', t, h2_pr, t, h2_po);
    grid on; ylabel('h_2 [cm]'); title(tytuly{ik});
    legend('nieliniowy', 'liniowy', sprintf('TS%d przed poprawą', L_wyb), ...
           sprintf('TS%d po poprawie', L_wyb), 'Location', 'best');

    subplot(2,2,ik+2);
    stairs(t, F1_dw{ik}); grid on; ylabel('F_{1in} [cm^3/s]'); xlabel('t [s]');
    ylim([0.4 1.6]*F1_0);
end
exportgraphics(f, fullfile(folder, 'zad2_model_przed_po_poprawie.png'), 'Resolution', 300); close(f);

%% Zapis funkcji przynależności wybranego modelu
save(plik_modelu, 'mf_wyb', 'mf_popr', 'h2_min', 'h2_max');
fprintf('\nZapisano funkcje przynależności do %s\n', plik_modelu);

%% Funkcje
% Wartości funkcji przynależności dla wektora h2 (wiersz = punkt, kolumna = zbiór)
function W = wartosci_mf(hh, mf)
    W = zeros(numel(hh), numel(mf.c));
    for k = 1:numel(hh)
        W(k,:) = zad2_przynaleznosc(hh(k), mf);
    end
end

% Modele lokalne w centrach zbiorów rozmytych
function lok = modele_lokalne(p, c)
    for i = numel(c):-1:1
        lok(i) = zad2_lokalny(p, c(i));
    end
end

% Model rozmyty TS: ważona suma pochodnych modeli lokalnych
function h2 = sym_ts(lok, mf, F1in, FD, Ts, d, F1_0, h1_0, h2_0)
    N = numel(F1in) - 1;
    h1 = h1_0; h2 = h2_0*ones(1, N+1);
    for k = 1:N
        if k > d, F1 = F1in(k-d); else, F1 = F1_0; end
        w = zad2_przynaleznosc(h2(k), mf);
        f1 = 0; f2 = 0;
        for i = 1:numel(lok)
            m = lok(i);
            f1 = f1 + w(i)*( m.d11*(h1 - m.h1) + m.b*(F1 + FD(k) - m.F) );
            f2 = f2 + w(i)*( m.d21*(h1 - m.h1) + m.d22*(h2(k) - m.h2) );
        end
        h1 = max(h1 + Ts*f1, 1e-6);
        h2(k+1) = max(h2(k) + Ts*f2, 0);
    end
end

% Model nieliniowy
function h2 = sym_nielin(p, F1in, FD, Ts, d, F1_0, h1_0, h2_0)
    N = numel(F1in) - 1;
    h1 = h1_0; h2 = h2_0*ones(1, N+1);
    for k = 1:N
        if k > d, F1 = F1in(k-d); else, F1 = F1_0; end
        F2 = p.a1*sqrt(h1); F3 = p.a2*sqrt(h2(k));
        h1 = max(h1 + Ts*(F1 + FD(k) - F2)/(2*p.C1*h1), 1e-6);
        h2(k+1) = max(h2(k) + Ts*(F2 - F3)/p.A2, 0);
    end
end
