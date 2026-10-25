% Znormalizowane wartości funkcji przynależności (suma = 1) dla wartości y = h2
function w = zad2_przynaleznosc(y, mf)
    c = mf.c;
    L = numel(c);
    if L == 1, w = 1; return; end
    y = min(max(y, c(1)), c(end));      % poza zakresem działa model skrajny
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
