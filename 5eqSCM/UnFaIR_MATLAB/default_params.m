function [gas_params , forcing_params , thermal_params] = default_params()

gas_params = [0.2173 1. 1. ; 0.2240, 0. 0. ; 0.2824 0. 0. ; 0.2763 0. 0. ;... % a1 : a4
              1000000. 9.15 116. ; 394.4 1. 1. ; 36.54 1. 1. ; 4.304 1. 1. ; ... % tau1 : tau4
              37.493303 8.54 67.2311 ; 0.01909 0. 0. ; 3.616153 -0.36 0. ; 0. 0.00031 -0.000906 ; ... % r0 : rA
              278.0 700.0 273.0 ]; % Preindustrial concentrations

gas_params(14,:) = 1./(5.148.*10^(18) ./ 1e18 .* [12.,16.,28.] ./ 28.97); % emission concentration conversion

forcing_params = [3.74/log(2) 0. 0. ; 0. 0. 0. ; 0 0.036 0.12]; % f1 : f3

thermal_params = [239. 4.1 ; 0.33 0.41]; % d ; q
          
end