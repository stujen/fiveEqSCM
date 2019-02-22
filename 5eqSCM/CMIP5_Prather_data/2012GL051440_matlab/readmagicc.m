function [year data unit] = readmagicc( rcpname, species, type)
%READMAGICC - Read emissions or concentrations for RCP scenarios.
%Emissions or concentration are from MAGICC6 model and include 1750-2500.
%
%
% Syntax:  [year,data,unit] = readmagicc( RCPname, species, type)
%
% Inputs:
%    RCPname - Name of RCP, as string. Options: RCP3PD, RCP45, RCP6, RCP85,
%              20THCENTURY
%    species - Name of species, as single string or cell array of strings. 
%              e.g. CH4, N2O, HFC134a
%    type    - 'E' for emissions, or 'C' for concentration output.
%
% Outputs:
%    year    - Column array of years
%    data    - Column array of emissions or concentration from MAGICC6
%              corresponding to elements of year
%    unit    - Cell array of string name of units for data
%
% Example 1: 
%    [year, emiss, unit] = readmagicc( 'RCP6', ['CH4';'N2O'], 'E' );
%    subplot(2,1,1)
%    plot( year, emiss(:,1) ); ylabel(strcat('CH4 emission, ',unit{1}))
%    subplot(2,1,2)
%    plot( year, emiss(:,2) ); ylabel(strcat('N2O emission, ',unit{2}))
%
% Example 2: 
%    [year, conc, unit] = readmagicc( 'RCP6', 'HFC134a', 'C' );
%    plot( year, conc ); ylabel(strcat('HFC134a concentration, ',unit{1}))
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Requires 'input/' subdirectory containing MAGICC6 emissions and mid-year
% concentrations in ASCII format, which can be downloaded from
% http://www.pik-potsdam.de/~mmalte/rcps/index.htm#Download
%

% Author: Christopher Holmes
% Department of Earth System Science
% University of California, Irvine
% email: cdholmes@uci.edu
% Website: http://www.ess.uci.edu/~cdholmes/
% November 2011; Last revision: 23-Feb-2012

%------------- BEGIN CODE --------------

% Choose emissions or concentrations
switch (lower(type))
    case ('e')
        % Filename
        file = 'input/%s_EMISSIONS.DAT';
        % Number of header lines
        headerlines=38;
        if strmatch(rcpname,'20THCENTURY') 
            headerlines=37;
        end
    case ('c')
        % Filename
        file = 'input/%s_MIDYEAR_CONCENTRATIONS.DAT';
        % Number of header lines
        headerlines=39;
end

% Change RCP names if necessary.
if strmatch(rcpname,'RCP60')
    filename = sprintf(file, 'RCP6');
elseif strmatch(rcpname,'RCP26')
    filename = sprintf(file, 'RCP3PD');
else
    filename = sprintf(file, rcpname);
end

% File is space delimited
DELIMITER = ' ';

% Import the file, ASCII format
newData = importdata(filename, DELIMITER, headerlines);

% Units and variable names for data in file
units = strread( char(newData.textdata(headerlines-1,1)), '%s', 'delimiter',' ');
names = newData.colheaders;

% Find the column for years
j=strmatch('YEARS',names);
year = newData.data(:,j);

% Convert to cell array of strings, if not already
species = cellstr( species );

% Number of requested species
D=length(species);

% Loop over requested species
for i=1:D

    % Find the column for requested species i
    j=strmatch(strtrim(char(species{i})),names);
    
    % Save the value and units for requested species i
    data(:,i) = newData.data(:,j);
    unit{i} = strtrim(char(units{j}));
  
end

end