functions {
#include cpl.stan
#include pgstat.stan
#include cpl_interval_fold.stan
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

  array[N_intervals, max(N_dets)] matrix[max_n_chan, max_n_echan] response;



  array[N_intervals, max(N_dets), max_n_chan] int mask;
  array[N_intervals,max(N_dets)] int N_channels_used;

  vector[N_intervals] dl;
  vector[N_intervals] z;


  // int N_gen_spectra;
  // vector[N_gen_spectra] model_energy;

  /* int N_correlation; */
  /* vector[N_correlation] model_correlation; */


}

transformed data {
  real x_r[0];
  int x_i[0];

  real kev2erg = 1.6021766e-9;
  real erg2kev = 6.24151e8;

  vector[N_intervals] dl2 = square(dl);
  int N_total_channels = 0;
  real emin = 10.;
  real emax = 1E4;

  array[N_intervals, max(N_dets)] vector[max_n_echan] ebounds_add;
  array[N_intervals, max(N_dets)] vector[max_n_echan] ebounds_half;

  array[N_intervals] int all_N;


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

  //vector<lower=-1.8, upper=1.>[N_intervals] alpha;
  vector<lower=-1.9, upper=1>[N_intervals] alpha;
  vector<lower=0, upper=4>[N_intervals] log_ec;
  //vector<lower=-5,upper=1>[N_intervals] log_K;

  //vector<lower=0, upper=5>[N_intervals] log_epeak;
  // vector<lower=0>[N_intervals] log_epeak;

  real log_energy_flux_mu_raw;
  real<lower=0> log_energy_flux_sigma;


  vector[N_intervals] log_energy_flux_raw; // raw energy flux norm

  

}

transformed parameters {

  vector[N_intervals] ec = pow(10, log_ec);
  vector[N_intervals] log_energy_flux;
  real log_energy_flux_mu;
  vector[N_intervals] energy_flux;

  vector[N_intervals] K;


  log_energy_flux_mu = log_energy_flux_mu_raw -7;
  
  log_energy_flux = (log_energy_flux_mu) + log_energy_flux_raw * log_energy_flux_sigma;
  energy_flux = pow(10, log_energy_flux);
  //vector[N_intervals] epeak;
  //vector[N_intervals] log_energy_flux;

  for (n in 1:N_intervals){

    array[3] real theta = {1., alpha[n], ec[n]};

    //epeak[n] = (2+alpha[n]) * pow(10, log_ec[n]);

    print(theta);
    K[n] = erg2kev * energy_flux[n]  * inv(integrate_1d(cpl_flux_integrand, 10., 1.e4, theta, x_r, x_i));
    //K[n] = erg2kev * energy_flux[n] * inv( ggrb_int_cpl(alpha[n], ec[n], 10., 1.e3) );

  }



}


model {

  int grainsize = 1;

  // alpha ~ normal(-1,.5);

  // log_epeak ~ normal(2.,1);

  //log_energy_flux ~ normal(-7, 1);

  log_energy_flux_mu_raw ~ std_normal();
  log_energy_flux_sigma ~ std_normal();
  
  alpha ~ normal(-1,.5);

  log_ec ~ normal(2.,1);

  // log_K ~ normal(-1, 1);

  // print(alpha);
  // print(log_ec);
  // print(log_K);

  target += reduce_sum(partial_log_like, all_N, grainsize,  alpha,  ec,  K,  observed_counts,  background_counts, background_errors,  mask, N_channels_used,exposure,  ebounds_lo,  ebounds_hi,  ebounds_add,  ebounds_half, response, idx_background_zero, idx_background_nonzero, N_bkg_zero, N_bkg_nonzero, N_dets,  N_chan,  N_echan,  max_n_chan,  emin,  emax) ;

}
