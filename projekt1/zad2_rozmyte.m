% Modele i regulatory rozmyte Takagi-Sugeno: funkcje przynależności dla zmiennej
% harmonogramującej h2, modele lokalne, symulacja modelu TS i projekt rozmytego DMC.

classdef zad2_rozmyte
    methods (Static)
        function mf = utworz_mf(L, typ, g, szer, h2_min, h2_max)
            x = linspace(0, 1, L);
            mf.L = L;                                % L - liczba modeli lokalnych
            mf.typ = typ;                            % typ - 'trojkatne' lub 'gauss'
            mf.g = g;                                % g - rozkład centrów: g = 1 podział równomierny, g > 1 centra gęściej przy małych h2
            mf.szer = szer;                          % szer - szerokość funkcji Gaussa względem odstępu między sąsiednimi centrami 
            mf.c = h2_min + (h2_max - h2_min)*x.^g;  % c - centra zbiorów (punkty linearyzacji modeli lokalnych)
        end

        function mf = mf_liniowy(h2_0)
            mf = struct('L', 1, 'typ', 'trojkatne', 'g', 1, 'szer', 0, 'c', h2_0);
        end

        function w = przynaleznosc(y, mf)
            c = mf.c;
            L = numel(c);
            if L == 1, w = 1; return; end
            y = min(max(y, c(1)), c(end));  % poza zakresem działa model skrajny
            switch mf.typ
                case 'trojkatne'
                    w = zeros(1, L);
                    for i = 1:L
                        if i > 1 && y >= c(i-1) && y <= c(i)
                            w(i) = (y - c(i-1))/(c(i) - c(i-1));
                        end
                        if i < L && y >= c(i) && y <= c(i+1)
                            w(i) = (c(i+1) - y)/(c(i+1) - c(i));
                        end
                    end
                case 'gauss'
                    % szerokość każdej funkcji zależy od odległości do sąsiednich centrów
                    dc  = diff(c);
                    odl = [dc(1), (dc(1:end-1) + dc(2:end))/2, dc(end)];
                    sig = mf.szer*odl;
                    w = exp(-(y - c).^2./(2*sig.^2));
                otherwise
                    error('Nieznany typ funkcji przynależności: %s', mf.typ);
            end
            w = w/sum(w);
        end

        function W = wartosci(hh, mf)
            W = zeros(numel(hh), numel(mf.c));
            for k = 1:numel(hh)
                W(k,:) = zad2_rozmyte.przynaleznosc(hh(k), mf);
            end
        end

        % Warianty kształtu do strojenia
        function war = warianty(L, g_lista, szer_lista, h2_min, h2_max)
            war = cell(1, numel(g_lista)*(1 + numel(szer_lista)));
            iw = 0;
            for g = g_lista
                iw = iw + 1;
                war{iw} = zad2_rozmyte.utworz_mf(L, 'trojkatne', g, 0, h2_min, h2_max);
                for sz = szer_lista
                    iw = iw + 1;
                    war{iw} = zad2_rozmyte.utworz_mf(L, 'gauss', g, sz, h2_min, h2_max);
                end
            end
        end

        % Krótki opis funkcji przynależności do tabel i tytułów wykresów
        function s = opis(mf)
            if strcmp(mf.typ, 'gauss')
                s = sprintf('%d zb., gauss, g = %g, szer = %g', numel(mf.c), mf.g, mf.szer);
            else
                s = sprintf('%d zb., trojkatne, g = %g', numel(mf.c), mf.g);
            end
        end

        % Modele lokalne zlinearyzowane w centrach zbiorów rozmytych
        function modele_lok = modele_lokalne(model, c)
            for i = numel(c):-1:1
                modele_lok(i) = zad1_obiekt.linearyzacja(model, c(i));
            end
        end

        % Symulacja modelu Takagi - Sugeko
        function h2 = sym_ts(modele_lok, mf, F1in, FD, Ts, d, F1_0, h1_0, h2_0)
            N = numel(F1in) - 1;
            h1 = h1_0; h2 = h2_0*ones(1, N+1);
            for k = 1:N
                if k > d, F1 = F1in(k-d); else, F1 = F1_0; end
                w = zad2_rozmyte.przynaleznosc(h2(k), mf);
                f1 = 0; f2 = 0;
                for i = 1:numel(modele_lok)
                    model_lokalny = modele_lok(i);
                    f1 = f1 + w(i)*( model_lokalny.d11*(h1 - model_lokalny.h1) + model_lokalny.b*(F1 + FD(k) - model_lokalny.F) );
                    f2 = f2 + w(i)*( model_lokalny.d21*(h1 - model_lokalny.h1) + model_lokalny.d22*(h2(k) - model_lokalny.h2) );
                end
                h1 = max(h1 + Ts*f1, 1e-6);
                h2(k+1) = max(h2(k) + Ts*f2, 0);
            end
        end

        % Rozmyty DMC
        function regulator = projekt_dmc(model, mf, nastawy, Ts, co, d)
            modele_lok = zad2_rozmyte.modele_lokalne(model, mf.c);
            regulator = zad1_dmc.projekt(modele_lok, nastawy, Ts, co, d);
            regulator.mf   = mf;
            regulator.wagi = @(y) zad2_rozmyte.przynaleznosc(y, mf);
        end

        % Wykres funkcji przynależności
        function wykres_mf(hh, mf_lista, tytuly, uklad, plik)
            if uklad(1) == 1
                f = zad1_wykresy.nowy(550*uklad(2), 400);
            else
                f = zad1_wykresy.nowy();
            end
            for i = 1:numel(mf_lista)
                subplot(uklad(1), uklad(2), i);
                plot(hh, zad2_rozmyte.wartosci(hh, mf_lista{i})); grid on;
                xlabel('h_2 [cm]'); ylabel('\mu'); title(tytuly{i});
            end
            zad1_wykresy.zapisz(f, plik);
        end
    end
end
