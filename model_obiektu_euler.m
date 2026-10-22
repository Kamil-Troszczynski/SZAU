%% Parametry
p.A2 = 430; p.C1 = 0.85; p.a1 = 18; p.a2 = 24; p.tau = 70;
F1_0 = 55; FD_0 = 10;
h1_0 = ((F1_0 + FD_0)/p.a1)^2;
h2_0 = ((F1_0 + FD_0)/p.a2)^2;

%% Dyskretyzacja
Ts = 0.5;                 % krok [s]
Tk = 3000;                % czas symulacji [s]
N  = round(Tk/Ts);
d  = round(p.tau/Ts);     % opóźnienie w próbkach

t    = (0:N)*Ts;
h1   = zeros(1,N+1);  h2 = zeros(1,N+1);
F1in = F1_0*ones(1,N+1);
h1(1) = h1_0;  h2(1) = h2_0;

%% Pętla symulacji
for k = 1:N
    F1in(k) = F1_0 + 5*(t(k) >= 100);
    FD      = FD_0 + 5*(t(k) >= 1500);

    if k > d
        F1 = F1in(k-d);
    else
        F1 = F1_0;
    end

    F2 = p.a1*sqrt(h1(k));
    F3 = p.a2*sqrt(h2(k));

    h1(k+1) = h1(k) + Ts*(F1 + FD - F2)/(2*p.C1*h1(k));
    h2(k+1) = h2(k) + Ts*(F2 - F3)/p.A2;

    h1(k+1) = max(h1(k+1), 1e-6);
    h2(k+1) = max(h2(k+1), 0);
end

%% Wykresy
subplot(2,1,1); plot(t, h1); grid on; ylabel('h_1 [cm]');
subplot(2,1,2); plot(t, h2); grid on; ylabel('h_2 [cm]'); xlabel('t [s]');