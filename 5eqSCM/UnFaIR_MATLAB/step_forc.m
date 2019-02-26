function RF = step_forc(C,PI_conc,f)

RF = f(1,:) .* log(C ./ PI_conc) + f(2,:) .* (C - PI_conc) + f(3,:) .* ( sqrt(C) - sqrt(PI_conc) );

end