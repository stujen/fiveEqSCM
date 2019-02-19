function input = readparams( species )

% Filename with parameters
filename = sprintf('input/inputparams%s.txt', species );

% Read parameters
fid = fopen( filename, 'r' );

line = textscan(fid,'%s',1, 'delimiter', '\n'); 

% Split string at spaces
subs = strread( char(line{1}), '%s', 'delimiter',' ');

% first entry is number of species in file
nS = str2num( char( subs(1) ) );

% initialize all variables with 0 uncertainty
Year   = 0;
Xair   = 0;
Mass        = zeros(nS,1);
cPD_glob    = zeros(nS,2);
dcdt_glob   = zeros(nS,2);
cPI_glob    = zeros(nS,2);
kOH_glob    = zeros(nS,2);
kCl_glob    = zeros(nS,2);
kStrat_glob = zeros(nS,2);
kSurf_glob  = zeros(nS,2);
anPI_glob   = [ones(nS,1), zeros(nS,1)];
aPI_glob    = [ones(nS,1), zeros(nS,1)];
a2100_glob  = [ones(nS,1), zeros(nS,1)];
aPIstrat_glob = [ones(nS,1), zeros(nS,1)];
a2100strat_glob = [ones(nS,1), zeros(nS,1)];
sOH_glob    = zeros(nS,2);
sStrat_glob = zeros(nS,2);
b_glob      = zeros(nS,2);
RFe_glob    = zeros(nS,2);
kMCF_glob   = zeros(nS,2);
kMCFstrat_glob = zeros(nS,2);
kMCFocean_glob = zeros(nS,2);
fillMCF_glob   = [ones(nS,1), zeros(nS,1)];
fill_glob = [ones(nS,1), zeros(nS,1)];
r272_glob   = zeros(nS,2);
r225_glob   = zeros(nS,2);

% species number
s=1;

% Read line with species name
line = textscan(fid,'%s',1, 'delimiter', '\n'); 

% Split string at spaces
subs = strread( char(line{1}), '%s', 'delimiter',' ');

% first entry is name of gas species
species = char( subs{1} );

% Skip headers
line = textscan(fid,'%s',1, 'delimiter', '\n'); 

while (~feof(fid))
               
    % Read one line
    line = textscan(fid,'%s',1, 'delimiter', '\n'); 
       
    % Split string at spaces
    subs = strread( char(line{1}), '%s', 'delimiter',' ');

    % Check for end of record for a single species
    if (strmatch('--',char(subs(1))))
        % Increment species counter
        s=s+1;
    
        if (s > nS)
            break
        end
        
        % Read line with species name
        line = textscan(fid,'%s',1, 'delimiter', '\n');
        
        % Split string at spaces
        subs = strread( char(line{1}), '%s', 'delimiter',' ');
        
        % first entry is name of gas species
        species = char( species, char(subs{1}) );       
        
        % skip 1 line, then read 2nd
        line = textscan(fid,'%s',1, 'delimiter', '\n');
        line = textscan(fid,'%s',1, 'delimiter', '\n');
        
        % Split string at spaces
        subs = strread( char(line{1}), '%s', 'delimiter',' ');
    
    end  
    
    % 1st entry is variable name
    name = char( subs(1) );
    
    % 2nd and 3rd entries are numerical values
    values = str2num( char( subs(2:3) ) )';
    
    % 3rd entry is units
    unit = char( subs(4) );
    
    % final entries are description
    desc = strcat( char( subs(5:end) ) );
    
    switch name
        case 'Year'
            Year = values(1);
        case 'Mass'
            Mass(s) = values(1);
        case 'Xair'
            Xair = values(1);
        case 'cPD'
            cPD_glob(s,:) = values;
        case 'dcdt'
            dcdt_glob(s,:) = values;
        case 'cPI'
            cPI_glob(s,:) = values;
        case 'kOH'
            kOH_glob(s,:) = values;
        case 'kCl'
            kCl_glob(s,:) = values;
        case 'kStrat'
            kStrat_glob(s,:) = values;
        case 'kSurf'
            kSurf_glob(s,:) = values;
        case 'anPI'
            anPI_glob(s,:) = values;
        case 'aPI'
            aPI_glob(s,:) = values;
        case 'a2100'
            a2100_glob(s,:) = values;
        case 'aPIstrat'
            aPIstrat_glob(s,:) = values;
        case 'a2100strat'
            a2100strat_glob(s,:) = values;
        case 'sOH'
            sOH_glob(s,:) = values;
        case 'sStrat'
            sStrat_glob(s,:) = values;
        case 'b'
            b_glob(s,:) = values;
        case 'RFe'
            RFe_glob(s,:) = values;
        case 'kMCF'
            kMCF_glob(s,:) = values;
        case 'kMCFstrat'
            kMCFstrat_glob(s,:) = values;
        case 'kMCFocean'
            kMCFocean_glob(s,:) = values;
        case 'fillMCF'
            fillMCF_glob(s,:) = values;
        case 'fill'
            fill_glob(s,:) = values;
        case 'r272'
            r272_glob(s,:) = values;
        case 'r225'
            r225_glob(s,:) = values;
    end
end

fclose(fid);

% Read emissions from file
[YErcp85_glob Ercp85_glob unit] = readmagicc('RCP85',species,'e');
[YErcp60_glob Ercp60_glob unit] = readmagicc('RCP60',species,'e');
[YErcp45_glob Ercp45_glob unit] = readmagicc('RCP45',species,'e');
[YErcp26_glob Ercp26_glob unit] = readmagicc('RCP26',species,'e');

% Read concentrations from file
[YCrcp85_glob Crcp85_glob unit] = readmagicc('RCP85',species,'c');
[YCrcp60_glob Crcp60_glob unit] = readmagicc('RCP60',species,'c');
[YCrcp45_glob Crcp45_glob unit] = readmagicc('RCP45',species,'c');
[YCrcp26_glob Crcp26_glob unit] = readmagicc('RCP26',species,'c');


input = struct( 'Year', Year, 'species', species, 'mass', Mass, 'Xair', Xair,...
       'cPD', @cPD, 'dcdt', @dcdt, 'cPI', @cPI, ...
       'kOH', @kOH, 'kCl', @kCl, 'kStrat', @kStrat, 'kSurf', @kSurf, ...
       'anPI', @anPI, 'aPI', @aPI, 'a2100', @a2100, 'aPIstrat', @aPIstrat, 'a2100strat', @a2100strat, ...
       'sOH', @sOH, 'sStrat', @sStrat, 'b', @b, 'RFe', @RFe, ...
       'kMCF', @kMCF, 'kMCFstrat', @kMCFstrat, 'kMCFocean', @kMCFocean, 'fillMCF', @fillMCF, ...
       'fill', @fill, 'r272', @r272, 'r225', @r225, ...
       'Ercp26', @Ercp26, 'Ercp45', @Ercp45, 'Ercp60', @Ercp60, 'Ercp85', @Ercp85, ...
       'Crcp26', @Crcp26, 'Crcp45', @Crcp45, 'Crcp60', @Crcp60, 'Crcp85', @Crcp85 );
   
%% Define functions 

    function c = cPD(varargin)
        % Present-day concentration, ppb
        
        %%global cPD_glob;
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        c = cPD_glob(:,1) + r .* cPD_glob(:,2);
    end

    function c = cPI(varargin)
        % Preindustrial concentration, ppb
        
        %%global cPI_glob;
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        c = cPI_glob(:,1) + r .* cPI_glob(:,2);
    end

    function dcdt = dcdt(varargin)
        % Present-day growth rate, ppb/a
        
        %%global dcdt_glob;
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        dcdt = dcdt_glob(:,1) + r .* dcdt_glob(:,2);
    end


    function k = kMCF(varargin)
        % Present-day MCF decay frequency, 1/a
        
        %%global kMCF_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        k = kMCF_glob(:,1) + r .* kMCF_glob(:,2);
    end

    function k = kMCFstrat(varargin)
        % Present-day MCF loss frequency due to stratosphere, 1/a
        
        %%global kMCFstrat_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        k = kMCFstrat_glob(:,1) + r .* kMCFstrat_glob(:,2);
    end

    function k = kMCFocean(varargin)
        % Present-day MCF loss frequency due to ocean uptake, 1/a
        
        %%global kMCFocean_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        k = kMCFocean_glob(:,1) + r .* kMCFocean_glob(:,2);
    end

    function f = fillMCF(varargin)
        % MCF fill factor, %%global/troposphere mean mixing ratios, unitless
        
        %%global fillMCF_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        f = fillMCF_glob(:,1) + r .* fillMCF_glob(:,2);
    end

    function f = fill(varargin)
        % species fill factor, %%global/troposphere mean mixing ratios, unitless
        
        %%global fill_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        f = fill_glob(:,1) + r .* fill_glob(:,2);
    end

    function kk = r272(varargin)
        % k(Species+OH)/k(MCF+OH) at 272 K (includes uncertainty in temperature)
        
        %%global r272_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        kk = r272_glob(:,1) + r .* r272_glob(:,2);
    end

    function kk = r225(varargin)
        % k(Species+OH)/k(MCF+OH) at 272 K (includes uncertainty in temperature)
        
        %%global r225_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        kk = r225_glob(:,1) + r .* r225_glob(:,2);
    end

    function k = kOH()
        % Present-day loss frequency due to tropospheric OH, 1/a
        % derived quantity
        
        k = r272() ./ fill() .* fillMCF() .* ( kMCF() - kMCFstrat() - kMCFocean() );
        
    end

    function k = kCl(varargin)
        % Present-day loss frequency due to tropospheric Cl, 1/a
        
        %%global kCl_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        k = kCl_glob(:,1) + r .* kCl_glob(:,2);
    end

    function k = kStrat(varargin)
        % Present-day loss frequency due to all stratospheric processes, 1/a
        
        %%global kStrat_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        k = kStrat_glob(:,1) + r .* kStrat_glob(:,2);
    end

    function k = kSurf(varargin)
        % Present-day loss frequency due to surface deposition, 1/a
        
        %%global kSurf_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        k = kSurf_glob(:,1) + r .* kSurf_glob(:,2);
    end

    function a = anPI(varargin)
        % Ratio of PI/PD natural emissions, unitless
        
        %%global anPI_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        a = anPI_glob(:,1) + r .* anPI_glob(:,2);
    end

    function a = aPI(varargin)
        % Ratio of PI/PD loss frequency due to tropospheric OH, unitless
        
        %%global aPI_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        a = aPI_glob(:,1) + r .* aPI_glob(:,2);
    end

    function a = a2100(varargin)
        % Ratio of y2100/PD loss frequency due to tropospheric OH, unitless
        
        %%global a2100_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        a = a2100_glob(:,1) + r .* a2100_glob(:,2);
    end

    function a = aPIstrat(varargin)
        % Ratio of PI/PD loss frequency due to stratosphere, unitless
        
        %%global aPIstrat_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        a = aPIstrat_glob(:,1) + r .* aPIstrat_glob(:,2);
    end


    function a = a2100strat(varargin)
        % Ratio of y2100/PD loss frequency due to stratosphere, unitless
        
        %%global a2100strat_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        a = a2100strat_glob(:,1) + r .* a2100strat_glob(:,2);
    end

    function s = sOH(varargin)
        % Present-day OH feedback, unitless
        % (dln(OH)/dln(C))=(dln(kOH)/dln(C))
        
        %%global sOH_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        s = sOH_glob(:,1) + r .* sOH_glob(:,2);
    end

    function s = sStrat(varargin)
        % Present-day stratosphere feedback, unitless
        % (dln(kStrat)/dln(C))
        
        %%global sStrat_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        s = sStrat_glob(:,1) + r .* sStrat_glob(:,2);
    end

    function b = b(varargin)
        % Conversion between tropospheric abundance and total burden, Tg/ppb
        
        %global b_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        b = b_glob(:,1) + r .* b_glob(:,2);
    end

    function e = RFe(varargin)
        % Radiative forcing efficiency, W/m2/ppb
        
        %global RFe_glob
        
        if (nargin == 1)
            r = varargin{1};
        else
            r = randn(nS,1);
        end
        
        e = RFe_glob(:,1) + r .* RFe_glob(:,2);
    end

    function E=Ercp85(t)
        % RCP85 emissions, Tg(CH4)/a
        
        E = interp1( YErcp85_glob, Ercp85_glob, t );
        
    end

    function E=Ercp60(t)
        % RCP60 emissions, Tg(CH4)/a
        
        E = interp1( YErcp60_glob, Ercp60_glob, t );
              
    end

    function E=Ercp26(t)
        % RCP26 emissions, Tg(CH4)/a
        E = interp1( YErcp26_glob, Ercp26_glob, t );
    end

    function E=Ercp45(t)
        % RCP45 emissions, Tg(CH4)/a
        
        E = interp1( YErcp45_glob, Ercp45_glob, t );
        
    end

    function c=Crcp45(t)
        % RCP45 abundance, ppb
        
        c = interp1( YCrcp45_glob, Crcp45_glob, t );
        
    end

    function c=Crcp26(t)
        % RCP26 abundance, ppb
        
        c = interp1( YCrcp26_glob, Crcp26_glob, t );
        
    end

    function c=Crcp60(t)
        % RCP60 abundance, ppb
        
        c = interp1( YCrcp60_glob, Crcp60_glob, t );
        
    end

    function c=Crcp85(t)
        % RCP85 abundance, ppb
        
        c = interp1( YCrcp85_glob, Crcp85_glob, t );
               
    end
   
end

