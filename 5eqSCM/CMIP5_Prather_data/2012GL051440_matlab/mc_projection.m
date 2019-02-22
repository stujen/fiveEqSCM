function mc_projection( species, RCPname, Nmc, doProj )

fprintf(1,'Monte carlo %s %s, %i iterations\n', species, RCPname, Nmc );

param=readparams(species);

% Emissions and Abundance from MAGICC
[tErcp_glob Ercp_glob unit] = readmagicc(RCPname,species,'E');
[tCrcp_glob Crcp_glob unit] = readmagicc(RCPname,species,'C');

% Initialize
Epd_anthro = zeros(Nmc,1); % present-day anthropogenic emissions, Tg/a
Epd_nat    = zeros(Nmc,1); % present-day natural emissions, Tg/a
Epd = zeros(Nmc,1); % present-day total emissions, Tg/a
Epi = zeros(Nmc,1); % preindustrial total emissions, Tg/a
kMCFtotal = zeros(Nmc,1);
kOHpd = zeros(Nmc,1);
kStratPD = zeros(Nmc,1); 
%cF      = cell(Nmc,1);
%kFoh    = cell(Nmc,1);
%kFstrat = cell(Nmc,1);
if (doProj) 
    NT = 10;  % Number of times in interval 2010-2100
    tt = linspace(param.Year,2100,NT); % times to plot
    cF     = zeros(Nmc,NT); % projected abundances, ppb
    Ercp   = zeros(Nmc,NT); % rescaled rcp emissions, Tg/a
    kFoh   = zeros(Nmc,NT); % projected loss frequency due to trop OH, 1/a
    kFstrat= zeros(Nmc,NT); % projected stratospheric loss frequency, 1/a
    kF     = zeros(Nmc,NT); % projected total loss frequency, 1/a
    RFrcp  = zeros(Nmc,NT); % projected anthropogenic radiative forcing, W/m2
end

for i=1:Nmc

    % fill factor
    fill = param.fill();
    
    % present-day total loss of MCF, 1/a
    kMCFtotal(i) = param.fillMCF() * ...
        ( param.kMCF() - param.kMCFstrat() - param.kMCFocean() );
    
    % present-day loss due to OH, 1/a
    kOHpd(i)  = kMCFtotal(i) * param.r272() / fill; 
    
    % present-day loss due to stratosphere, 1/a
    kStratPD(i) = param.kStrat(); 
    
    % present-day loss due to other processes, 1/a
    kOther = param.kCl() + param.kSurf();
    
    % present-day loss frequency, 1/a
    kPD = kOHpd(i) + kStratPD(i) + kOther;
    
    % preindustrial loss frequency, 1/a
    kPI = param.aPI() * kOHpd(i) + param.aPIstrat() * kStratPD(i) + kOther;
       
    % Conversion factor ppb -> Tg
    ppb2Tg = param.Xair * param.mass * fill; 
    
    % present-day abundance, ppb
    cPD = param.cPD();
    
    % preindustrial abundance, ppb
    cPI = param.cPI();
     
    % present-day burden, Tg
    Bpd = ppb2Tg * cPD;
    
    % present-day emissions, Tg/a
    Epd(i) = kPD * Bpd + ppb2Tg * param.dcdt();
    
    % preindustrial burden, Tg
    Bpi = ppb2Tg * cPI;
    
    % preindustrial emissions, Tg/a
    Epi(i) = kPI * Bpi;
    
    % present-day natural emissions, Tg/a
    Epd_nat(i) = Epi(i) * param.anPI();
    
    % present-day anthropogenic emissions, Tg/a
    Epd_anthro(i) = Epd(i) - Epd_nat(i);

    if (Epd_anthro(i)<0)
        x=1;
    end
    
    if (doProj)
        % Future RCP emissions, scaled to match constraints on 2010, Tg/a
        E = @(t) Ercp0(t) * Epd_anthro(i) / Ercp0(param.Year);
        
        % Future change in loss frequency due to OH, unitless
        a = param.a2100();
        aF = @(t) interp1( [param.Year 2100], [1 a], t);
        
        % Future change in loss frequency due to stratosphere, unitless
        a = param.a2100strat();
        aFstrat = @(t) interp1( [param.Year 2100], [1 a], t);
        
        % present-day tropospheric OH feedback factor, dln(OH)/dln(C) = dln(kOH)/dln(C)
        sPD = param.sOH();
        
        % present-day stratospheric feedback factor,  dln(kStrat)/dln(C)
        sStratPD = param.sStrat();
        
        % Future loss frequency, 1/a
        %kF = @(t,c) ( kOHpd(i) * exp(sPD * log(c./cPD)) + kStratPD(i) + kOther ); % Includes trop OH feedback
        kF1 = @(t,c) kOHpd(i)    * exp( sPD * log(c./cPD) ) .* aF(t);
        kF2 = @(t,c) kStratPD(i) * exp( sStratPD * log(c./cPD) ) .* aFstrat(t);
        kFtot = @(t,c) kF1(t,c) + kF2(t,c) + kOther; % Includes trop and strat feedback and uncertainty in baseline OH
        
        % Define integrand
        dydt = @(t,c) ( E(t) + Epi(i) ) / ppb2Tg - kFtot(t,c) * c;
        
        % Integrate future burden, ppb
        sol = ode23( dydt, [param.Year,2100], cPD );
        
        % Future abundance at decadal intervals, ppb
 %       cF    = deval( sol, tt );
        
        % Save future abundance, ppb
        cF(i,:) = deval( sol, tt );
    
        % Save rescaled future RCP emissions, Tg/a
        Ercp(i,:) = E(tt);
    
        kFoh(i,:) = kF1(tt,cF(i,:));
        kFstrat(i,:) = kF2(tt,cF(i,:));
        kF(i,:) = kFtot(tt,cF(i,:));
        
        % Save future loss by trop OH, 1/a
%         kOHrcp(i,:) = kFoh(tt,Crcp(i,:));
%         
%         % Save future total loss, 1/a
%         krcp(i,:)   = kF(tt,Crcp(i,:));
        
        % Save future RF, W/m2
        RFrcp(i,:) = ( cF(i,:) - cPI ) * param.RFe();
    end
end

% Save monte carlo iterations
fname = sprintf( 'output/output_%s_%s_%d.mat', species, RCPname, Nmc ); 
save( fname, 'Epi', 'Epd', 'Epd_anthro', 'Epd_nat', 'param' );
if (doProj)
    save( fname, 'tt', 'Ercp', 'cF', 'RFrcp', 'kF', 'kFoh', 'kFstrat', 'kMCFtotal','-append');
%    save( fname, 'tt', 'Ercp', 'Crcp', 'RFrcp', 'krcp', 'kOHrcp', 'cF', 'kFoh', 'kFstrat', '-append');
end


    function E=Ercp0(t)
        % RCP emissions, Tg/a
        
        E = interp1( tErcp_glob, Ercp_glob, t );
              
    end


    function c=Crcp0(t)
        % RCP abundance, ppb
        
        c = interp1( tCrcp_glob, Crcp_glob, t );
               
    end

end