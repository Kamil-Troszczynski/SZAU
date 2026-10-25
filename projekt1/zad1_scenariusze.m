% Przebiegi wejść używane w eksperymentach

classdef zad1_scenariusze
    properties (Constant)
        T_F1   = [300 1200];    % czasy włączenia sterowania F1in
        M_F1   = [1.5 0.5];     % mnożniki M_<wielkość>, dla których testujemy skuteczność modeli
        T_YZAD = [300 1200];    % czsy ustawienia wartości zadanej h2
        M_YZAD = [1.5 0.5];
        T_FD   = [700 2000];    % czas włączenia zakłócenia FD
        M_FD   = [0.5 1.5];
    end

    methods (Static)
        function v = przebieg(t, x0, czasy, mnozniki)
            v = x0*ones(size(t));
            for i = 1:numel(czasy)
                v(t >= czasy(i)) = mnozniki(i)*x0;
            end
        end

        function F1in = przebieg_F1(t, F1_0)
            F1in = zad1_scenariusze.przebieg(t, F1_0, zad1_scenariusze.T_F1, zad1_scenariusze.M_F1);
        end

        function FD = przebieg_FD(t, FD_0)
            FD = zad1_scenariusze.przebieg(t, FD_0, zad1_scenariusze.T_FD, zad1_scenariusze.M_FD);
        end

        function yzad = przebieg_yzad(t, h2_0)
            yzad = zad1_scenariusze.przebieg(t, h2_0, zad1_scenariusze.T_YZAD, zad1_scenariusze.M_YZAD);
        end

        function F1_dw = wzrost_spadek(t, F1_0)
            F1_dw = {zad1_scenariusze.przebieg(t, F1_0, [100 1500], [1.25 1.5]), zad1_scenariusze.przebieg(t, F1_0, [100 1500], [0.75 0.5])};
        end

        function [F1_s, FD_s, opis] = skoki_pojedyncze(t, F1_0, FD_0, skoki_F1, skoki_FD, t_skok)
            nF1  = numel(skoki_F1);
            nS   = nF1 + numel(skoki_FD);
            F1_s = cell(1, nS);
            FD_s = cell(1, nS);
            opis = cell(1, nS);
            for is = 1:nS
                if is <= nF1
                    F1_s{is} = zad1_scenariusze.przebieg(t, F1_0, t_skok, 1 + skoki_F1(is));
                    FD_s{is} = FD_0*ones(size(t));
                    opis{is} = sprintf('F1 %+g%%', 100*skoki_F1(is));
                else
                    sk = skoki_FD(is - nF1);
                    F1_s{is} = F1_0*ones(size(t));
                    FD_s{is} = zad1_scenariusze.przebieg(t, FD_0, t_skok, 1 + sk);
                    opis{is} = sprintf('FD %+g%%', 100*sk);
                end
            end
        end
    end
end
