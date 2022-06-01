real ggrb_int_cpl(real alpha, real ec, real emin, real emax) {
  real i1 = gamma_q(2 + alpha, emin / ec) * tgamma(2 + alpha);
  real i2 = gamma_q(2 + alpha, emax / ec) * tgamma(2 + alpha);
  return -square(ec) * (i2-i1);
}



vector cpl(vector energy, real norm, real ec, real alpha) {

  real piv = 100.;

  // vector[num_elements(energy)] log_v = alpha * log(energy/piv) - energy/ec;

  // return norm * exp(log_v);

    
    
    return norm * pow(energy/piv, alpha) .* exp(-energy/ec);
  
}


real cpl_indi(real energy, real K, real alpha, real ec) {

  real piv = 100.;

  return K * pow(energy/piv, alpha) * exp(-energy/ec);
  
 
}

real cpl_flux_integrand(real x, real xc, real[] theta, real[] x_r, int[] x_i) {

  real out = x * cpl_indi(x, theta[1], theta[2], theta[3]);
  
  return out;


}



vector differential_flux( vector energy, real norm, real ec, real alpha) {

  return cpl(energy, norm, ec, alpha);
	  
}



vector integral_flux(vector ebounds_lo, vector ebounds_hi, vector ebounds_add, vector ebounds_half, real norm, real ec, real alpha) {

  return (ebounds_add
	  .* (differential_flux(ebounds_lo, norm, ec, alpha)
	      + 4 * differential_flux(ebounds_half, norm, ec, alpha)
	      + differential_flux(ebounds_hi, norm, ec, alpha)));
      
}
