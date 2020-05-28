real ggrb_int_pl(real alpha, real beta, real ec, real emin, real emax) {

  real pre = pow(alpha - beta, alpha - beta) * exp(beta - alpha) / pow(ec, beta);
  if (beta !=-2) {
    return pre/(2.+beta) * (pow(emax, 2+ beta) - pow(emin, 2+ beta));
  }
  else {
    return pre * log(emax/emin);
  }
}

  
real ggrb_int_cpl(real alpha, real ec, real emin, real emax) {
  real i1 = gamma_q(2 + alpha, emin / ec) * tgamma(2 + alpha);
  real i2 = gamma_q(2 + alpha, emax / ec) * tgamma(2 + alpha);
  return -square(ec) * (i2-i1);
}


real [] band_precalculation(real flux, real alpha, real  beta, real epeak, real emin, real emax) {

  real ec;
  real esplit;
  real intflux;
  real norm;
  real pre;
  real erg2keV = 6.24151e8;
  real alpha_minus_beta = alpha - beta;
    
  if (alpha !=-2.) {
    ec = epeak / (2. + alpha); 
  }
  else {
    ec = epeak/.0001; 
  }

  esplit = (alpha_minus_beta) * ec;
    
  if ((emin<= esplit)  &&  (esplit <=emax)) {
    intflux = (ggrb_int_cpl(alpha, ec, emin, esplit) + ggrb_int_pl(alpha, beta, ec, esplit, emax));
  }
    
  else if (esplit < emin) {
    intflux = ggrb_int_pl(alpha, beta, ec, esplit, emax);  
  }
    
  else {
    print(esplit);
  }
 
  norm = flux * erg2keV / intflux;
  pre = pow(alpha_minus_beta, alpha_minus_beta) * exp(-alpha_minus_beta);
  return {norm, ec, esplit, pre};
}

vector differential_flux( vector energy, real norm , real ec, real esplit, real alpha, real beta, real pre) {

  int N = num_elements(energy); 
  vector[N] out;
  vector[N] ratio;
       
  ratio = energy/ec;
    
  for (n in 1:N) {
      
    if (energy[n] < esplit) {
      out[n] = pow(ratio[n], alpha) * exp(-ratio[n]);
    }   
    else {
      out[n] = pre * pow(ratio[n], beta);
    }      
      
  }
    
  return norm * out;
	  
}

    
vector integral_flux(vector ebounds_lo, vector ebounds_hi, vector ebounds_add, vector ebounds_half,  real norm, real ec, real esplit, real alpha, real  beta, real pre) {

  return (ebounds_add
	  .* (differential_flux(ebounds_lo, norm, ec, esplit, alpha, beta, pre)
	      + 4 * differential_flux(ebounds_half, norm, ec, esplit, alpha, beta, pre)
	      + differential_flux(ebounds_hi, norm, ec, esplit, alpha, beta, pre  )));
      
}
