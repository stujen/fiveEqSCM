function g1 = g_1( a , tau , h )

if nargin < 3
    h = 100;
end

g1 = sum(a .* tau .* ( 1. - ( 1. + h./tau ) .* exp(-h./tau) ), 1 );

end