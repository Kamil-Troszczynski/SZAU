% Krótki opis funkcji przynależności do tabel i legend
function s = zad2_opis_mf(mf)
    if strcmp(mf.typ, 'gauss')
        s = sprintf('%d zb., gauss, g = %g, szer = %g', numel(mf.c), mf.g, mf.szer);
    else
        s = sprintf('%d zb., trojkatne, g = %g', numel(mf.c), mf.g);
    end
end
