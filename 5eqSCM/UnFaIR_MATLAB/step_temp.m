function [S,T] = step_temp(S,F,q,d)

S = q .* F .* ( 1 - exp(-1./d) ) + S .* exp(-1./d);

T = sum(S);

end