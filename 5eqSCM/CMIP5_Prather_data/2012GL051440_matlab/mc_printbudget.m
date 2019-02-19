function mc_printbudget( species, RCPname, N )

% Load monte carlo output
fname = sprintf( 'output/output_%s_%s_%d.mat', species, RCPname, N ); 
load( fname );

nS = size( param.species );
nS = nS(1);

for s=1:nS
    
    % species name
    species = param.species(s,:);
    
    if (mean(Epd(:,s)) < 1)
       unitA='ppt';
       unitE='Gg/a';
       unitRF='mW/m2';
       Epd(:,s) = Epd(:,s) * 1e3;
       Epi(:,s) = Epi(:,s) * 1e3;
       Epd_anthro(:,s) = Epd_anthro(:,s) * 1e3;
       Epd_nat(:,s) = Epd_nat(:,s) * 1e3;
       cF(:,:,s) = cF(:,:,s) * 1e3;
       Ercp(:,:,s) = Ercp(:,:,s) * 1e3;
       RFrcp(:,:,s) = RFrcp(:,:,s) * 1e3;
    else        
       unitA='ppb';
       unitE='Tg/a';
       unitRF='W/m2';
    end
    
    % Basic output
    disp(' ')
    disp(' ')
    disp('=================')
    fprintf(1,'%s %s\n',RCPname,species)
    disp('=================')
    fprintf(1,'%12s %8.2f +/- %8.2f\n', 'PI emiss', mean(Epi(:,s)), std(Epi(:,s)));
    fprintf(1,'%12s %8.2f +/- %8.2f\n', 'Y2010 emiss', mean(Epd(:,s)), std(Epd(:,s)));
    fprintf(1,'%12s %8.2f +/- %8.2f\n', 'Y2010 anthro', mean(Epd_anthro(:,s)), std(Epd_anthro(:,s)) );
    fprintf(1,'%12s %8.2f +/- %8.2f\n', 'Y2010 nature', mean(Epd_nat(:,s)), std(Epd_nat(:,s)) );
    disp(' ')
    
    if (exist('tt','var'))
        % Display projected abundance
        fprintf(1,'%4s  %10s(%s) %12s(%s) %10s(%s) %17s %21s\n', ...
            'Year','Emission',unitE,'Abundance',unitA, 'RF',unitRF, 'T_OH(a)', 'T_total(a)' )
        for t=1:length(tt)
            fprintf(1,'%4i %8.1f +/- %5.1f %8.1f +/- %5.1f %8.3f +/- %5.3f %8.2f +/- %5.2f %8.2f +/- %5.2f\n', ...
                tt(t), ...
                mean(Ercp(:,t,s)), std(Ercp(:,t,s)), ... % mean +/- sd emissions
                mean(cF(:,t,s)), std(cF(:,t,s)), ... % mean +/- sd abundance
                mean(RFrcp(:,t,s)), std(RFrcp(:,t,s)), ... % mean +/- sd RF
                1./mean(kFoh(:,t,s)), std(kFoh(:,t,s))/mean(kFoh(:,t,s))^2, ... % mean +/- sd lifetime due to trop OH
                1./mean(kF(:,t,s)),   std(kF(:,t,s))  /mean(kF(:,t,s))^2 ) % mean +/- sd total lifetime
        end
    end
    
end

if (nS > 1)
    disp(' ')
    disp('All species in file')
    fprintf(1,'%4s  %11s(%s)\n', 'Year','RF', unitRF )
    for t=1:length(tt)
        fprintf(1,'%4i %8.3f +/- %5.3f\n', ...
            tt(t), ...
            mean( sum(RFrcp(:,t,:),3) ), std( sum(RFrcp(:,t,:),3) ) )% mean +/- sd RF
    end
end

end