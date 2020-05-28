real partial_log_like(int [] n_slice, int start, int end, vector alpha, vector log_epeak, vector log_energy_flux, vector[,] observed_counts, vector[,] background_counts, vector[,] background_errors, int [,,] mask, int [,] N_channels_used, real [,] exposure, vector[,] ebounds_lo, vector[,] ebounds_hi, vector[,] ebounds_add, vector[,] ebounds_half, matrix[,] response, int [,,] idx_background_zero, int [,,] idx_background_nonzero, int [,] N_bkg_zero, int [,] N_bkg_nonzero, int [] N_dets, int [,] N_chan, int [,] N_echan, int max_n_chan, real emin, real emax ) {


  real log_like = 0;
  int slice_length = num_elements(n_slice);
  real pre_calc[2];

  vector[max_n_chan] expected_model_counts[max(N_dets)];

  for (i in 1:slice_length) {

    int n = n_slice[i];

    // norm, ec, epslit, pre
    pre_calc[:] = band_precalculation(10^log_energy_flux[n], alpha[n], 10^log_epeak[n], emin, emax);



    for (m in 1:N_dets[n]) {


      expected_model_counts[m, : N_chan[n,m]] = ((to_row_vector(integral_flux(ebounds_lo[n, m, :N_echan[n, m]],
                                                                              ebounds_hi[n, m, :N_echan[n, m]],
                                                                              ebounds_add[n, m, :N_echan[n, m]],
                                                                              ebounds_half[n, m, :N_echan[n, m]],
                                                                              pre_calc[1],
                                                                              pre_calc[2],
                                                                              alpha[n]
                                                                              )) * response[n, m,:N_echan[n,m],:N_chan[n,m]]))' * exposure[n,m];



      log_like += sum(pgstat(observed_counts[n, m, mask[n,m,:N_channels_used[n,m]]],
                             background_counts[n, m, mask[n,m,:N_channels_used[n,m]]],
                             background_errors[n, m, mask[n,m,:N_channels_used[n,m]]],
                             expected_model_counts[m, mask[n,m,:N_channels_used[n,m]]],
                             idx_background_zero[n,m, :N_bkg_zero[n,m]],
                             idx_background_nonzero[n,m, :N_bkg_nonzero[n,m]]  ));

    }

  }

  return log_like;
}
