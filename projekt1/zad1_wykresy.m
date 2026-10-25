% Wspólne wykresy: tworzenie i zapis figur, przebiegi regulacji, porównanie modeli

classdef zad1_wykresy
    methods (Static)
        function f = nowy(szer, wys)
            if nargin < 2
                f = figure('Color', 'w', 'Visible', 'off');
            else
                f = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 szer wys]);
            end
        end

        function zapisz(f, plik)
            exportgraphics(f, plik, 'Resolution', 300);
            close(f);
        end

        % Przebiegi regulacji
        function regulacja(t, ip, yzad, FD, h2, u, nazwy_h2, nazwy_u, plik)
            f = zad1_wykresy.nowy(1000, 700);
            subplot(2,1,1); hold on; box on;
            for i = 1:numel(h2)
                stairs(t(ip), h2{i}(ip));
            end
            stairs(t(ip), yzad(ip), 'k--'); hold off;
            grid on; ylabel('h_2 [cm]');
            legend([nazwy_h2, {'wartość zadana'}], 'Location', 'best');

            subplot(2,1,2); hold on; box on;
            for i = 1:numel(u)
                stairs(t(ip), u{i}(ip));
            end
            if isempty(FD)
                ylabel('F_{1in} [cm^3/s]');
            else
                stairs(t(ip), FD(ip), 'k--');
                ylabel('F [cm^3/s]');
                legend([nazwy_u, {'zakłócenie F_D'}], 'Location', 'best');
            end
            hold off; grid on; xlabel('t [s]');
            zad1_wykresy.zapisz(f, plik);
        end

        % Wykres 2x2: h2 (góra) i F1in (dół) dla wzrostu (lewa) i spadku (prawa) F1in.
        % h2{ik} - cell z przebiegami dla kierunku ik, style{i} - styl linii ('' = domyślny)
        function wzrost_spadek(t, F1_dw, F1_0, h2, nazwy, style, plik)
            tytuly = {'Wzrost F_{1in}', 'Spadek F_{1in}'};
            f = zad1_wykresy.nowy(1200, 700);
            for ik = 1:2
                subplot(2,2,ik); hold on; box on;
                for i = 1:numel(h2{ik})
                    if isempty(style{i})
                        plot(t, h2{ik}{i});
                    else
                        plot(t, h2{ik}{i}, style{i});
                    end
                end
                hold off; grid on; ylabel('h_2 [cm]'); title(tytuly{ik});
                legend(nazwy, 'Location', 'best');

                subplot(2,2,ik+2);
                stairs(t, F1_dw{ik}); grid on; ylabel('F_{1in} [cm^3/s]'); xlabel('t [s]');
                ylim([0.4 1.6]*F1_0);
            end
            zad1_wykresy.zapisz(f, plik);
        end
    end
end
