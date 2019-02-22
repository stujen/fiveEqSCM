function mc_script
% script to run monte carlo simulations

%number of random realizations
N = 100000;

% do projections? if false, then program will just calculate PI and PD budgets
doProj = true;


mc_projection( 'N2O', 'RCP26', N, doProj )
mc_projection( 'N2O', 'RCP45', N, doProj )
mc_projection( 'N2O', 'RCP60', N, doProj )
mc_projection( 'N2O', 'RCP85', N, doProj )

mc_printbudget( 'N2O', 'RCP26', N )
mc_printbudget( 'N2O', 'RCP45', N )
mc_printbudget( 'N2O', 'RCP60', N )
mc_printbudget( 'N2O', 'RCP85', N )

mc_projection( 'CH4', 'RCP26', N, doProj )
mc_projection( 'CH4', 'RCP45', N, doProj )
mc_projection( 'CH4', 'RCP60', N, doProj )
mc_projection( 'CH4', 'RCP85', N, doProj )

mc_printbudget( 'CH4', 'RCP26', N )
mc_printbudget( 'CH4', 'RCP45', N )
mc_printbudget( 'CH4', 'RCP60', N )
mc_printbudget( 'CH4', 'RCP85', N )

if (doProj)
    mc_projection_hfc( 'HFC134a', 'RCP26', N )
    mc_projection_hfc( 'HFC134a', 'RCP45', N )
    mc_projection_hfc( 'HFC134a', 'RCP60', N )
    mc_projection_hfc( 'HFC134a', 'RCP85', N )

    mc_printbudget( 'HFC134a', 'RCP26', N )
    mc_printbudget( 'HFC134a', 'RCP45', N )
    mc_printbudget( 'HFC134a', 'RCP60', N )
    mc_printbudget( 'HFC134a', 'RCP85', N )

%     % All other HFCs
%     mc_projection_hfc( 'HFCall', 'RCP26', N )
%     mc_projection_hfc( 'HFCall', 'RCP45', N )
%     mc_projection_hfc( 'HFCall', 'RCP60', N )
%     mc_projection_hfc( 'HFCall', 'RCP85', N )
%     mc_printbudget( 'HFCall', 'RCP26', N )
%     mc_printbudget( 'HFCall', 'RCP45', N )
%     mc_printbudget( 'HFCall', 'RCP60', N )
%     mc_printbudget( 'HFCall', 'RCP85', N )

    mc_plotprojection( 'CH4', N )
    mc_plotprojection( 'N2O', N )
    mc_plotprojection( 'HFC134a', N )
end

end
