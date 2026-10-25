% Model lokalny obiektu: linearyzacja w punkcie równowagi o poziomie h2i
%   dh1/dt = d11*(h1 - m.h1) + b*(F1 + FD - m.F)
%   dh2/dt = d21*(h1 - m.h1) + d22*(h2 - m.h2)
function m = zad2_lokalny(p, h2i)
    m.h2  = h2i;
    m.F   = p.a2*sqrt(h2i);             % całkowity dopływ F1 + FD w równowadze
    m.h1  = (m.F/p.a1)^2;
    m.d11 = -p.a1/(4*p.C1*m.h1^(3/2));
    m.b   =  1/(2*p.C1*m.h1);
    m.d21 =  p.a1/(2*p.A2*sqrt(m.h1));
    m.d22 = -p.a2/(2*p.A2*sqrt(h2i));
end
