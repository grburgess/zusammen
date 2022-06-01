functions {
#include pgstat.stan
#include band_grb.stan
}

data {

  int<lower=1> N_intervals;
  int max_n_echan;
  int max_n_chan;

  array[N_intervals] int<lower=0> N_dets; // number of detectors poer data type
  array[N_intervals, max(N_dets)] int<lower=0> N_chan; // number of channels in each detector
  array[N_intervals,  max(N_dets)] int<lower=0> N_echan; // number of energy side channels in each detector

  array[N_intervals] int grb_id;
  int N_grbs;
  
  array[N_intervals, max(N_dets)] vector[max_n_echan] ebounds_hi;
  array[N_intervals, max(N_dets)] vector[max_n_echan] ebounds_lo;
  


  array[N_intervals, max(N_dets)] vector[max_n_chan] observed_counts;
  array[N_intervals, max(N_dets)] vector[max_n_chan] background_counts;
  array[N_intervals, max(N_dets)] vector[max_n_chan] background_errors;

  array[N_intervals, max(N_dets), max_n_chan] int idx_background_zero;
  array[N_intervals, max(N_dets), max_n_chan] int idx_background_nonzero;
  array[N_intervals,max(N_dets)] int N_bkg_zero;
  array[N_intervals,max(N_dets)] int N_bkg_nonzero;
  
  array[N_intervals, max(N_dets)] real exposure;

  array[N_intervals, max(N_dets)] matrix[max_n_echan, max_n_chan] response;


  
  array[N_intervals, max(N_dets), max_n_chan] int mask;
  array[N_intervals,max(N_dets)] int N_channels_used;

  vector[N_intervals] dl;
  vector[N_intervals] z;


  int N_gen_spectra;
  vector[N_gen_spectra] model_energy;

  
  
}

transformed data {
  vector[N_intervals] dl2 = square(dl); 
  

  int N_total_channels = 0;
  real emin = 10.;
  real emax = 1E6;
  
  

  
  array[N_intervals, max(N_dets)] vector[max_n_echan] ebounds_add;
  array[N_intervals, max(N_dets)] vector[max_n_echan] ebounds_half;
  
  for (n in 1:N_intervals) {

    for (m in 1:N_dets[n]) {

      ebounds_half[n, m, :N_echan[n, m]] = 0.5*(ebounds_hi[n, m, :N_echan[n, m]]+ebounds_lo[n, m, :N_echan[n, m]]);
      ebounds_add[n, m, :N_echan[n, m]] = (ebounds_hi[n, m, :N_echan[n, m]] - ebounds_lo[n, m, :N_echan[n, m]])/6.0;
      N_total_channels += N_channels_used[n,m];
    }
    
    
  }

      
  /* vector[N_intervals] emin; */
  /* vector[N_intervals] emax; */

  /* emin = 10. ./ (1+z); */
  /* emax = 1.E5 ./ (1+z); */
  

  
}



parameters {

  vector<lower=-1.8, upper=1.>[N_intervals] alpha;
  vector<lower=-6., upper=-2.>[N_intervals] beta;
  vector<lower=0, upper=4>[N_intervals] log_epeak;
  vector[N_intervals] log_energy_flux;
  
}

transformed parameters {

  array[N_intervals, max(N_dets)] vector[max_n_chan] expected_model_counts;
  real pre_calc[N_intervals, 4]; 
  
  // compute the folded counts
  
    for (n in 1:N_intervals) {

      // norm, ec, epslit, pre 
      pre_calc[n, :] = band_precalculation(10^log_energy_flux[n], alpha[n], beta[n], 10^log_epeak[n], emin, emax);
      
      for (m in 1:N_dets[n]) {
	
	expected_model_counts[n,m,:N_chan[n,m]] = ((to_row_vector(integral_flux(ebounds_lo[n, m, :N_echan[n, m]],
										ebounds_hi[n, m, :N_echan[n, m]],
										ebounds_add[n, m, :N_echan[n, m]],
										ebounds_half[n, m, :N_echan[n, m]],
										pre_calc[n,1],
										pre_calc[n,2],
										pre_calc[n,3],
										alpha[n],
										beta[n],
										pre_calc[n,4])) * response[n, m,:N_echan[n,m],:N_chan[n,m]]) * exposure[n,m])';
	
      }
    }
  
  
}


model {

  vector[N_total_channels] log_like;
  int pos = 1;
  
  alpha ~ normal(-1,.5);
  beta ~ normal(-3,1);
  log_epeak ~ normal(2.,1.);
  log_energy_flux ~ normal(-6,2);
  
  for (n in 1:N_intervals) {

    for (m in 1:N_dets[n]) {
      
      
      log_like[pos : pos + N_channels_used[n,m] - 1] = pgstat(observed_counts[n, m, mask[n,m,:N_channels_used[n,m]]],
							  background_counts[n, m, mask[n,m,:N_channels_used[n,m]]],
							  background_errors[n, m, mask[n,m,:N_channels_used[n,m]]],
							  expected_model_counts[n,m, mask[n,m,:N_channels_used[n,m]]],
							  idx_background_zero[n,m, :N_bkg_zero[n,m]],
							  idx_background_nonzero[n,m, :N_bkg_nonzero[n,m]]);
      
      pos += N_channels_used[n,m];
      
    }
    
  }

  target += sum(log_like);
  
  
  
}


generated quantities {

  vector[N_gen_spectra] vfv_spectra[N_intervals];
  vector[max_n_chan] count_ppc[N_intervals, max(N_dets)];
  vector[max_n_chan] source_ppc[N_intervals, max(N_dets)];
  

  for (n in 1:N_intervals) {

    vfv_spectra[n] =square(model_energy) .* differential_flux(model_energy, pre_calc[n, 1], pre_calc[n, 2], pre_calc[n, 3], alpha[n], beta[n], pre_calc[n, 4]);


      
 
    

    for (m in 1:N_dets[n]) {

        vector[N_channels_used[n,m]] ppc_background = background_model(observed_counts[n, m, mask[n,m,:N_channels_used[n,m]]],
								       background_counts[n, m, mask[n,m,:N_channels_used[n,m]]],
								       background_errors[n, m, mask[n,m,:N_channels_used[n,m]]],
								       expected_model_counts[n, m, mask[n,m,:N_channels_used[n,m]]]);
	
	vector[N_channels_used[n,m]] rate = ppc_background + expected_model_counts[n, m, mask[n,m,:N_channels_used[n,m]]] ;
	vector[N_channels_used[n,m]] source_rate = expected_model_counts[n, m, mask[n,m,:N_channels_used[n,m]]] ;

	
	for (i in 1:N_channels_used[n,m]) {

	  if (rate[i]>2^30) {
	    

	    count_ppc[n,m,i] = 0;
	    
	  }
	  
	  else {

	    
	    count_ppc[n,m,i] = poisson_rng( rate[i] );
	    
	  }



	  if (source_rate[i]>2^30) {
	    
	    
	    source_ppc[n,m,i] = 0;
	    
	  }
	  
	  else {
	    
	    
	    source_ppc[n,m,i] = poisson_rng( source_rate[i] );
	    
	  }
	  
	  
	}
	
    }
    
    
    
  
  }    
  
  
  
}

  
