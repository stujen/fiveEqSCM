function alpha = alpha_val(G,G_A,T,a,tau,r,h,iirf100_max)

if nargin < 8
    iirf100_max = 97.;
    if nargin < 7
        h = 100.;
    end
end

iirf100 = r(1,:) + r(2,:) .* (G-G_A) + r(3,:) .* T + r(4,:) .* G_A; % r0 + rC * (G-G_A) + rT * T + rA * G_A

iirf100 = abs(iirf100); % take the absolute value

for i = 1:length(iirf100)
    if iirf100(i) > iirf100_max
        iirf100(i) = iirf100_max; % if iirf100 is larger than maximum allowed value, cap it
    end
end

alpha = g_0(a , tau , h) .* sinh( iirf100 ./ g_1(a , tau , h) ); % calculate alpha from iirf100 & g0, g1

end
    