function mc_projection_hfc( species, RCPname, Nmc )
% This requires that CH4 projections have already been done in order to
% calculate future OH

fprintf(1,'Monte carlo %s %s, %i iterations\n', species, RCPname, Nmc );

% Load monte carlo output for CH4
fname = sprintf( 'output/output_%s_%s_%d.mat', 'CH4', RCPname, Nmc ); 
ch4=load( fname );

%species = 'HFC134a';
param = readparams( species );
nHFC = size(param.species,1);

[yr em unit] = readmagicc( RCPname, param.species, 'e' );
E_HFC = @(t) interp1( yr, em, t)/1d3;

tt=ch4.tt;
NT = length( tt );
kMCFtotal = ch4.kMCFtotal;

Ercp   = zeros(Nmc,NT,nHFC);
cF     = zeros(Nmc,NT,nHFC);
kF     = zeros(Nmc,NT,nHFC);
kFoh   = zeros(Nmc,NT,nHFC);
kFstrat= zeros(Nmc,NT,nHFC);
RFrcp  = zeros(Nmc,NT,nHFC);

for i=1:Nmc
    
    % fill factor, unitless
    fill = param.fill();
    
    % present-day loss due to OH, 1/a
    kOHpd = ch4.kMCFtotal(i) * param.r272() ./ fill; % for HFCs
   
    % Conversion factor ppb -> Tg
    ppb2Tg = param.Xair * param.mass .* fill; 

    % Future loss due to OH, 1/a
    k1 = @(t) kOHpd/ ch4.kFoh(i,1) * interp1( ch4.tt, ch4.kFoh(i,:), t );
    
    % Future loss due to stratosphere, 1/a
    k2 = @(t) param.r225() * interp1( ch4.tt, ch4.kFstrat(i,:), t );
    
%    kOH = @(t) kOHpd/ ch4.kOHpd(i) * ch4.kFoh{i}( t, ch4.cF{i}(t) );
%    kStrat = @(t) kStratPD / ch4.kStratPD(i) * ch4.kFstrat{i}( t, ch4.cF{i}(t) );
    %kF = @(t) kOH(t) + kStrat(t);
    kFtotal = @(t) k1(t) + k2(t);

    % Define integrand for HFCs
    dydt =  @(t,c) E_HFC(t)' ./ ppb2Tg - kFtotal(t) .* c ;

    % Integrate future burden, ppb
    sol = ode23( dydt, [2010,2100], param.cPD() );    

    % Future CH4 abundance at decadal intervals, ppb
    %cF(i,:,:) = deval( sol, tt );
    
    % Save future RF , W/m2    
    RFe_hfc = param.RFe();
    for t=1:NT
        cF(i,t,:) = deval( sol, tt(t) );
        % Save future loss by trop OH, 1/a
        kFoh(i,t,:) = k1(tt(t));
        kFstrat(i,t,:)=k1(tt(t));
        kF(i,t,:)  = kFtotal(tt(t));
        RFrcp(i,t,:) = squeeze(cF(i,t,:)) .* RFe_hfc;
        Ercp(i,t,:) = E_HFC(tt(t));
    end
end

Epi = zeros(1,nHFC);
Epd = Ercp(:,1,:);
Epd_nat = zeros(1,nHFC);
Epd_anthro = Epd;

% Save monte carlo iterations
fname = sprintf( 'output/output_%s_%s_%d.mat', species, RCPname, Nmc ); 
save( fname, 'Epi', 'Epd', 'Epd_anthro', 'Epd_nat', 'param' );
save( fname, 'tt', 'Ercp', 'cF', 'RFrcp', 'kF', 'kFoh', 'kFstrat', 'kMCFtotal','-append');

end