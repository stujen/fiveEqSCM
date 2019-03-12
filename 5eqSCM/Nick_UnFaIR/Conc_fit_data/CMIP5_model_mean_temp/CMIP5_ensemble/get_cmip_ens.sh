for i in {01..38};
do
  echo $i
  wget "https://climexp.knmi.nl/data/icmip5_tas_Amon_mod_rcp85_0-360E_-90-90N_n_0"$i"_mean1_anom_30.dat"
done
