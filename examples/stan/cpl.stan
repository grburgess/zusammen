  
real ggrb_int_cpl(real alpha, real ec, real emin, real emax) {
  real i1 = gamma_q(2 + alpha, emin / ec) * tgamma(2 + alpha);
  real i2 = gamma_q(2 + alpha, emax / ec) * tgamma(2 + alpha);
  return -square(ec) * (i2-i1);
}


real [] band_precalculation(real flux, real alpha, real epeak, real emin, real emax) {

  real ec;
  real intflux;
  real norm;
  
  real erg2keV = 6.24151e8;
  
    
  if (alpha !=-2.) {
    ec = epeak / (2. + alpha); 
  }
  else {
    ec = epeak/.0001; 
  }

    
  intflux = ggrb_int_cpl(alpha, ec, emin, emax);
 
    
  
  norm = flux * erg2keV / intflux;

  return {norm, ec};
}

vector differential_flux( vector energy, real norm , real ec, real alpha) {

  int N = num_elements(energy); 
  vector[N] out;
  vector[N] ratio;
       
  ratio = energy/ec;
    
  for (n in 1:N) {
      

    out[n] = norm * pow(ratio[n], alpha) * exp(-ratio[n]);
  }
    
  return out;
	  
}



vector integral_flux(vector ebounds_lo, vector ebounds_hi, vector ebounds_add, vector ebounds_half, real norm, real ec, real alpha) {

  return (ebounds_add
	  .* (differential_flux(ebounds_lo, norm, ec, alpha)
	      + 4 * differential_flux(ebounds_half, norm, ec, alpha)
	      + differential_flux(ebounds_hi, norm, ec, alpha)));
      
}
