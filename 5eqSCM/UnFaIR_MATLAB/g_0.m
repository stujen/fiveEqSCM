function g0 = g_0( a , tau , h )

if nargin < 3
    h = 100;
end

g0 = ( sinh( sum( a .* tau .* ( 1. - exp(-h./tau) ) , 1) ./ g_1(a,tau,h) ) ).^(-1.);

end