% Model obiektu - parametry, linearyzacja i symulacja

classdef zad1_obiekt
    methods (Static)
        % Parametry obiektu i punkt pracy
        function [model, F1_0, FD_0, h1_0, h2_0] = parametry()
            model.A2 = 430; model.C1 = 0.85; model.a1 = 18; model.a2 = 24; model.tau = 70;
            F1_0 = 55; FD_0 = 10;
            h1_0 = ((F1_0 + FD_0)/model.a1)^2;
            h2_0 = ((F1_0 + FD_0)/model.a2)^2;
        end

        % Model liniowy
        function model_zlinearyzowany = linearyzacja(model, h2r)
            model_zlinearyzowany.h2 = h2r;
            model_zlinearyzowany.F = model.a2*sqrt(h2r);
            model_zlinearyzowany.h1 = (model_zlinearyzowany.F/model.a1)^2;
            model_zlinearyzowany.d11 = -model.a1/(4*model.C1*model_zlinearyzowany.h1^(3/2));
            model_zlinearyzowany.b = 1/(2*model.C1*model_zlinearyzowany.h1);
            model_zlinearyzowany.d21 = model.a1/(2*model.A2*sqrt(model_zlinearyzowany.h1));
            model_zlinearyzowany.d22 = -model.a2/(2*model.A2*sqrt(h2r));
        end

        % Wzmocnienie statyczne modelu liniowego
        function K = wzmocnienie(model_zlinearyzowany)
            K = (model_zlinearyzowany.b/(-model_zlinearyzowany.d11))*(model_zlinearyzowany.d21/(-model_zlinearyzowany.d22));
        end

        % Jeden krok metody Eulera dla obiektu nieliniowego
        function [h1, h2] = krok(model, h1, h2, F1, FD, Ts)
            F2 = model.a1*sqrt(h1);
            F3 = model.a2*sqrt(h2);
            h1n = h1 + Ts*(F1 + FD - F2)/(2*model.C1*h1);
            h2n = h2 + Ts*(F2 - F3)/model.A2;
            h1 = max(h1n, 1e-6);
            h2 = max(h2n, 0);
        end

        % Symulacja obiektu nieliniowego z wykorzystaniem metody Eulera.
        % F1in - sterowanie (przed opóźnieniem d próbek), FD - zakłócenie
        function [h2, h1] = sym_nieliniowy(model, F1in, FD, Ts, d, F1_0, h1_0, h2_0)
            N = numel(F1in) - 1;
            h1 = h1_0*ones(1, N+1);  
            h2 = h2_0*ones(1, N+1);
            for k = 1:N
                if k > d
                    F1 = F1in(k-d); 
                else
                    F1 = F1_0;
                end
                [h1(k+1), h2(k+1)] = zad1_obiekt.krok(model, h1(k), h2(k), F1, FD(k), Ts);
            end
        end

        % Symulacja modelu liniowego
        function [h2, h1] = sym_liniowy(model_zlinearyzowany, F1in, FD, Ts, d, F1_0, h1_0, h2_0)
            N = numel(F1in) - 1;
            h1 = h1_0*ones(1, N+1);  h2 = h2_0*ones(1, N+1);
            for k = 1:N
                if k > d
                    F1 = F1in(k-d); 
                else 
                    F1 = F1_0; 
                end
                h1(k+1) = max(h1(k) + Ts*( model_zlinearyzowany.d11*(h1(k) - model_zlinearyzowany.h1) + model_zlinearyzowany.b*(F1 + FD(k) - model_zlinearyzowany.F) ), 1e-6);
                h2(k+1) = max(h2(k) + Ts*( model_zlinearyzowany.d21*(h1(k) - model_zlinearyzowany.h1) + model_zlinearyzowany.d22*(h2(k) - model_zlinearyzowany.h2) ), 0);
            end
        end

        % Prawe strony równań różniczkowych dla ode45
        function dh = ode(t, h, model, F1, FD)
            h1 = max(h(1), 1e-6);
            h2 = max(h(2), 0);
            F2 = model.a1*sqrt(h1);
            F3 = model.a2*sqrt(h2);
            dh = [ (F1(t) + FD(t) - F2) / (2*model.C1*h1) ; (F2 - F3) / model.A2 ];
        end
    end
end
