function [C,RF,T,alpha] = UnFaIR(emissions , F_ext , gas_params, forcing_params, thermal_params)

if nargin < 3
    [gas_params , forcing_params , thermal_params] = default_params(); % assume default params if none given
    if nargin < 2
        F_ext = 0.; % if not given, assume zero external forcing
    end
end

a = gas_params(1:4,:);
tau = gas_params(5:8,:);
r = gas_params(9:12,:);
PI_conc = gas_params(13,:);
emis2conc = gas_params(14,:);

f = forcing_params;

d = thermal_params(1,:);
q = thermal_params(2,:);

G = cumsum(emissions);
C = zeros(size(emissions));
RF = zeros(size(emissions));
T = zeros(length(emissions),1);
alpha = zeros(size(emissions));

if size(F_ext,1) ~= 1
    if size(emissions,1) ~= size(F_ext,1)
        disp("External forcing length is not the same as emissions length!")
        return
    end
else
    F_ext = ones(length(emissions),1) .* F_ext;
end

alpha(1,:) = alpha_val(0.,0.,0.,a,tau,r);
[C(1,:),R,G_A] = step_conc(zeros(4,3) , emissions(1,:) , alpha(1,:) , a , tau , PI_conc , emis2conc);
RF(1,:) = step_forc(C(1,:) , PI_conc , f);
[S,T(1)] = step_temp(zeros(1,2),sum(RF(1,:))+F_ext(1,:),q,d);

for t = 2:(length(emissions))
    alpha(t,:) = alpha_val(G(t-1,:),G_A,T(t-1),a,tau,r);
    [C(t,:),R,G_A] = step_conc(R , emissions(t,:) , alpha(t,:) , a , tau , PI_conc , emis2conc);
    RF(t,:) = step_forc(C(t,:) , PI_conc , f);
    [S,T(t)] = step_temp(S,sum(RF(t,:))+F_ext(t,:),q,d);
end

RF = horzcat(RF,F_ext);
RF = horzcat(RF,sum(RF,2));
    
end