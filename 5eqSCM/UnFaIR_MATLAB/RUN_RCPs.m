RCP = input ('choose RCP (3PD,45,6,85)');

ems = csvread(strcat('/home/nleach/Documents/UnFaIR/5eqSCM/RCP_data/RCP',RCP,'_EMISSIONS.csv'), 37);
emissions = zeros(length(ems),3);
emissions(:,1) = sum(ems(:,2:3),2);
emissions(:,2:3) = ems(:,4:5);

forc = csvread(strcat('/home/nleach/Documents/UnFaIR/5eqSCM/RCP_data/RCP',RCP,'_MIDYEAR_RADFORCING.csv'), 59);
otherforc = forc(:,2) - forc(:,8);

[C,RF,T,alpha] = UnFaIR(emissions,otherforc);

plot_var = input( 'choose variable to plot (C, RF, T, alpha)' );

if or(strcmp(plot_var,'C'),strcmp(plot_var,'alpha'))
    plot_gas = input( 'choose gas to plot (CO2, CH4, N2O)' );
end

if strcmp(plot_var,'C')
    if strcmp(plot_gas , 'CO2')
        plot(C(:,1))
    elseif strcmp(plot_gas , 'CH4')
        plot(C(:,2))
    elseif strcmp(plot_gas , 'N2O')
        plot(C(:,3))
    end
    
elseif strcmp(plot_var,'alpha')
    if strcmp(plot_gas , 'CO2')
        plot(alpha(:,1))
    elseif strcmp(plot_gas , 'CH4')
        plot(alpha(:,2))
    elseif strcmp(plot_gas , 'N2O')
        plot(alpha(:,3))
    end

elseif strcmp(plot_var,'T')
    plot(T)
    
elseif strcmp(plot_var,'RF')
    plot(RF)
end 
    