functions {
#include cpl.stan
#include pgstat.stan
#include cpl_interval_fold.stan
}

data {

  int<lower=1> N_intervals;
  int max_n_echan;
  int max_n_chan;

  int<lower=0> N_dets[N_intervals]; // number of detectors poer data type
  int<lower=0> N_chan[N_intervals, max(N_dets)]; // number of channels in each detector
  int<lower=0> N_echan[N_intervals,  max(N_dets)]; // number of energy side channels in each detector

  int grb_id[N_intervals];
  int N_grbs;

  vector[max_n_echan] ebounds_hi[N_intervals, max(N_dets)];
  vector[max_n_echan] ebounds_lo[N_intervals, max(N_dets)];



  vector[max_n_chan] observed_counts[N_intervals, max(N_dets)];
  vector[max_n_chan] background_counts[N_intervals, max(N_dets)];
  vector[max_n_chan] background_errors[N_intervals, max(N_dets)];

  int idx_background_zero[N_intervals, max(N_dets), max_n_chan];
  int idx_background_nonzero[N_intervals, max(N_dets), max_n_chan];
  int N_bkg_zero[N_intervals,max(N_dets)];
  int N_bkg_nonzero[N_intervals,max(N_dets)];

  real exposure[N_intervals, max(N_dets)];

  matrix[max_n_echan, max_n_chan] response[N_intervals, max(N_dets)];



  int mask[N_intervals, max(N_dets), max_n_chan];
  int N_channels_used[N_intervals,max(N_dets)];

  vector[N_intervals] dl;
  vector[N_intervals] z;


  int N_gen_spectra;
  vector[N_gen_spectra] model_energy;

  /* int N_correlation; */
  /* vector[N_correlation] model_correlation; */


}

transformed data {
  vector[N_intervals] dl2 = square(dl);
  int N_total_channels = 0;
  real emin = 10.;
  real emax = 1E4;
  vector[max_n_echan] ebounds_add[N_intervals, max(N_dets)];
  vector[max_n_echan] ebounds_half[N_intervals, max(N_dets)];

  int all_N[N_intervals];


  // precalculation of energy bounds

  for (n in 1:N_intervals) {

    all_N[n] = n;
    
    for (m in 1:N_dets[n]) {
      ebounds_half[n, m, :N_echan[n, m]] = 0.5*(ebounds_hi[n, m, :N_echan[n, m]]+ebounds_lo[n, m, :N_echan[n, m]]);
      ebounds_add[n, m, :N_echan[n, m]] = (ebounds_hi[n, m, :N_echan[n, m]] - ebounds_lo[n, m, :N_echan[n, m]])/6.0;
      N_total_channels += N_channels_used[n,m];
    }


  }




}



parameters {

  vector<lower=-1.8, upper=1.>[N_intervals] alpha;

  vector<lower=0, upper=5>[N_intervals] log_epeak;
  vector[N_intervals] log_energy_flux;





}

transformed parameters {

}


model {

  //vector[N_total_channels] log_like;
  //int pos = 1;
  int grainsize = 1;


  alpha ~ normal(-1,.5);

  log_epeak ~ normal(2.,1);

  log_energy_flux ~ normal(-7, 1);


  target += reduce_sum(partial_log_like, all_N, grainsize,  alpha,  log_epeak,  log_energy_flux,  observed_counts,  background_counts, background_errors,  mask, N_channels_used, exposure,  ebounds_lo,  ebounds_hi,  ebounds_add,  ebounds_half,  response,   idx_background_zero,   idx_background_nonzero,  N_bkg_zero, N_bkg_nonzero, N_dets,  N_chan,  N_echan,  max_n_chan,  emin,  emax) ;

}


generated quantities {

  /* vector[N_gen_spectra] vfv_spectra[N_intervals]; */
  // vector[max_n_chan] count_ppc[N_intervals, max(N_dets)];
  /* //  vector[max_n_chan] source_ppc[N_intervals, max(N_dets)]; */
  /* /\* vector[N_correlation] correlations[N_grbs]; *\/ */

  /* for (n in 1:N_intervals) { */
  /*   vfv_spectra[n] =square(model_energy) .* differential_flux(model_energy, pre_calc[n, 1], pre_calc[n, 2], alpha[n]); */

  /* for (m in 1:N_dets[n]) { */

  /*   /\* vector[N_channels_used[n,m]] ppc_background = background_model(observed_counts[n, m, mask[n,m,:N_channels_used[n,m]]], *\/ */
  /*   /\*                                                             background_counts[n, m, mask[n,m,:N_channels_used[n,m]]], *\/ */
  /*    /\*                                                            background_errors[n, m, mask[n,m,:N_channels_used[n,m]]], *\/ */
  /*    /\*                                                            expected_model_counts[n, m, mask[n,m,:N_channels_used[n,m]]]); *\/ */

  /*   //       vector[N_channels_used[n,m]] rate = ppc_background + expected_model_counts[n, m, mask[n,m,:N_channels_used[n,m]]] ; */
  /*   vector[N_channels_used[n,m]] source_rate = expected_model_counts[n, m, mask[n,m,:N_channels_used[n,m]]]; */

  /*   for (i in 1:N_channels_used[n,m]) { */


  /*    /\* if (rate[i]<=0) { *\/ */
  /*    /\*   print(pre_calc[n,:]); *\/ */
  /*    /\*   print(ppc_background); *\/         */
  /*    /\* } *\/                  */
  /*    /\* if (rate[i]>2^30) { *\/ */
  /*    /\*   count_ppc[n,m,i] = 0; *\/      */
  /*    /\* } *\/          */
  /*    /\* else { *\/       */
  /*    /\*   count_ppc[n,m,i] = poisson_rng( rate[i] ); *\/         */
  /*    /\* } *\/          */
  /*    if (source_rate[i]>2^30) { */
  /*      source_ppc[n,m,i] = 0; */
  /*    } */

  /*    else { */
  /*      source_ppc[n,m,i] = poisson_rng( source_rate[i] ); */
  /*    } */
  /*   } */
  /* } */

  // }


  /* for (n in 1:N_grbs) { */

  /*   correlations[n] = delta[n] + gamma[n] * (model_correlation + log_zp1_grb[n] - 2 ) - log_dl2_grb[n]; */

  /* } */

}
