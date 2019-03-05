pro UnFaIR, emissions, F_ext, gas_params=gas_params, forcing_params=forcing_params, thermal_params=thermal_params

  IF ~KEYWORD_SET(gas_params) THEN BEGIN
    gas_params = default_gas_params()
  ENDIF
  
  IF ~KEYWORD_SET(forcing_params) THEN BEGIN
    forcing_params = default_forcing_params()
  ENDIF
  
  IF ~KEYWORD_SET(thermal_params) THEN BEGIN
    thermal_params = default_thermal_params()
  ENDIF

  ; Set the parameter values

  a = gas_params(*,0:3)
  tau = gas_params(*,4:7)
  r = gas_params(*,8:11)
  PI_conc = gas_params(*,12)
  emis2conc = gas_params(*,13)
  
  f = forcing_params
  
  d = thermal_params(*,0)
  q = thermal_params(*,1)
  
  ; Create the variable arrays
  
  G = TOTAL(emissions,2,/cumulative)
  ems_shape = size(emissions,/dimensions)
  C = MAKE_ARRAY(ems_shape)
  RF = MAKE_ARRAY(ems_shape)
  T = MAKE_ARRAY(ems_shape(1))
  alpha = MAKE_ARRAY(ems_shape)
  
  ; initialize the first timestep
  
  G_A = MAKE_ARRAY(ems_shape(0))
  R_pools = MAKE_ARRAY(size(a,/dimensions))
  
  alpha(*,0) = alpha_val(0,0,0, a, tau, r)
  
  concstep_result = step_conc(R_pools , emissions(*,0) , alpha(*,0) , a , tau , PI_conc , emis2conc)
  C(*,0) = concstep_result(*,0)
  R_pools = concstep_result(*,1:-2)
  G_A = concstep_result(*,-1)
  
  RF(*,0) = step_forc(C(*,0), PI_conc, f)
  
  S = MAKE_ARRAY(2)
  
  tempstep_result = step_temp(S, TOTAL(RF) + F_ext(0), q , d)
  T(0) = tempstep_result(0)
  S = tempstep_result(1:-1)
  
  for time = 1 , ems_shape(1) - 1 DO begin
    
    alpha(*,time) = alpha_val(G(*,time-1), G_A, T(time-1), a, tau, r)

    concstep_result = step_conc(R_pools , emissions(*,time) , alpha(*,time) , a , tau , PI_conc , emis2conc)
    C(*,time) = concstep_result(*,0)
    R_pools = concstep_result(*,1:-2)
    G_A = concstep_result(*,-1)

    RF(*,time) = step_forc(C(*,time), PI_conc, f)

    tempstep_result = step_temp(S, TOTAL(RF(*,time)) + F_ext(time), q , d)
    T(time) = tempstep_result(0)
    S = tempstep_result(1:-1)
    
  endfor

  save, filename='Model_output.sav' , C,RF,T,alpha,emissions,R_pools,S  

end

pro RUN_RCPS, NONE=none

  RCP = ''
  READ, RCP, PROMPT='Choose RCP: RCP3PD, RCP45, RCP6, RCP85'
  emissions = GET_RCP_EMS(RCP)
  forcing = GET_RCP_FORC(RCP,'TOTAL_INCLVOLCANIC_RF') - GET_RCP_FORC(RCP,'CO2CH4N2O_RF')
  
  UnFaIR, emissions, forcing

end

function GET_RCP_EMS, RCP

  RCP_filename = '/home/nleach/Documents/UnFaIR/5eqSCM/RCP_data/'+RCP+'_EMISSIONS.csv'
  RCP_ems = Read_csv(RCP_filename,N_TABLE_HEADER=36,HEADER=RCP_gas_names)
  return, TRANSPOSE([[RCP_ems.FIELD02 + RCP_ems.FIELD03] , [RCP_ems.FIELD04] , [RCP_ems.FIELD05]])

end

function GET_RCP_FORC, RCP, FORC_NAME

  RCP_filename = '/home/nleach/Documents/UnFaIR/5eqSCM/RCP_data/'+RCP+'_MIDYEAR_RADFORCING.csv'
  RCP_forc = Read_csv(RCP_filename,N_TABLE_HEADER=58,HEADER=RCP_forc_names)
  
  FORC_INDEX = WHERE(STRMATCH(RCP_FORC_NAMES, FORC_NAME) EQ 1)
  
  return, RCP_FORC.(FORC_INDEX)

end

function step_conc, R_pools, E, alpha, a, tau, PI_conc, emis2conc

  shape = size(a,/dimensions)
  E = rebin(E,shape)
  emis2conc = rebin(emis2conc,shape)
  alpha = rebin(emis2conc,shape)

  R_pools = E * emis2conc * a * alpha * tau * ( 1. - exp( -1./(alpha*tau) ) ) + R_pools * exp( -1./(alpha*tau) )
  
  C = PI_conc + TOTAL(R_pools,2)
  
  G_A = (C - PI_CONC) / emis2conc
  
  return, [[C],[R_pools],[G_A]]

end

function step_forc, C, PI_conc, f
  
  RF = f(*,0) * alog(C/PI_conc) + f(*,1) * (C - PI_conc) + f(*,2) * (sqrt(C) - sqrt(PI_conc))
  
  return, RF

end

function step_temp, S, RF, q, d

  S = q * RF * ( 1. - exp(-1./d) ) + S * exp(-1./d)
  
  T = TOTAL(S)
  
  return, [T,S]

end

function default_gas_params,NONE=none

  return, [[0.2173, 1., 1. ],[ 0.2240, 0., 0. ],[ 0.2824, 0., 0. ],[ 0.2763, 0., 0.], $ ; a1 : a4
  [1000000., 9.15, 116. ],[ 394.4, 1., 1. ],[ 36.54, 1., 1. ],[ 4.304, 1., 1. ], $ ; tau1 : tau4
  [37.493303, 8.54, 67.2311],[ 0.01909, 0., 0.],[ 3.616153, -0.36, 0.],[ 0., 0.00031, -0.000906], $ ; r0 : rA
  [278.0, 700.0, 273.0], $ ; Preindustrial concentration
  [0.4690, 0.3517, 0.2010]] ; emission to concentration conversions

end

function default_forcing_params,NONE=none

  return, [[3.74/alog(2), 0., 0. ],[ 0., 0., 0. ],[ 0., 0.036, 0.12]] ; f1 : f3

end

function default_thermal_params,NONE=none

  return, [[239., 4.1 ],[ 0.33, 0.41 ]] ;  d , q

end

pro save_default_parameter_state, gas_params=gas_params, forcing_params=forcing_params, thermal_params=thermal_params

  IF ~KEYWORD_SET(gas_params) THEN BEGIN
    gas_params = default_gas_params()
  ENDIF
  
  IF ~KEYWORD_SET(forcing_params) THEN BEGIN
    forcing_params = default_forcing_params()
  ENDIF
  
  IF ~KEYWORD_SET(thermal_params) THEN BEGIN
    thermal_params = default_thermal_params()
  ENDIF

  a = gas_params(*,0:3)
  tau = gas_params(*,4:7)
  r = gas_params(*,8:11)
  PI_conc = gas_params(*,12)
  emis2conc = gas_params(*,13)

  f = forcing_params

  d = thermal_params(*,0)
  q = thermal_params(*,1)
  
  save, /variables , filename='UnFaIR_default_parameter_state.sav'

end

function g_1, a , tau , h = h

; This calculates the g_1 coefficients for each gas 

  IF ~KEYWORD_SET(h) THEN BEGIN
    h = 100
  ENDIF

  return, TOTAL( a * tau * ( 1. - ( 1.  + h/tau) * exp(-h / tau) ) , 2)

END


function g_0, a , tau , h = h

; this calculates the g_0 coefficients for each gas

  IF ~KEYWORD_SET(h) THEN BEGIN
    h = 100
  ENDIF
  
  return, (sinh( TOTAL( a * tau * (1. - exp(-h / tau)) , 2) / g_1(a,tau,h=h) ) )^(-1.)
  
END


function alpha_val, G,G_A,T,a,tau,r,h=h,iirf100_max=iirf100_max

; This calculates the fractional timescale adjustment based upon the model state

  IF ~KEYWORD_SET(h) THEN BEGIN
    h = 100
  ENDIF
  
  IF ~KEYWORD_SET(iirf100_max) THEN BEGIN
    iirf100_max = 97.
  ENDIF

  iirf100 = r[*,0] + r[*,1] * (G-G_A) + r[*,2] * T + r[*,3] * G_A
  
  iirf100 = abs(iirf100)
  
  iirf100 = iirf100 < iirf100_max
  
  RETURN, g_0(a,tau,h=h) * sinh(iirf100 / g_1(a,tau,h=h))

END




