function [C, R, G_A] = step_conc(R,E,alpha,a,tau,PI_conc,emis2conc)

R = E .* emis2conc .* a .* alpha .* tau .* ( 1 - exp( -1 ./ (alpha.*tau) ) ) + R .* exp( -1 ./ (alpha.*tau) );

C = PI_conc + sum(R,1);

G_A = (C - PI_conc) ./ emis2conc;

end