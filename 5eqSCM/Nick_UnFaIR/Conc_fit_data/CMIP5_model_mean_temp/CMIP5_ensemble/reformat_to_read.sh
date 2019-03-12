	for i in {01..38};
do
sed -n '/1861  /,$p' "icmip5_tas_Amon_mod_rcp85_0-360E_-90-90N_n_0"$i"_mean1_anom_30.dat" > "icmip5_tas_Amon_mod_rcp85_0-360E_-90-90N_n_0"$i"_mean1_anom_30_reformat.dat"
done
