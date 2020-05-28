vector background_model(vector observed_counts, vector background_counts, vector background_error, vector expected_model_counts) {
  
  int N = num_elements(expected_model_counts);
  
  vector[N] MB = background_counts + expected_model_counts;
  vector[N] s2 = square(background_error);
  
  vector[N] b = 0.5 * (sqrt(square(MB) - 2 * s2 .* (MB - 2 * observed_counts) + square(s2))
		       + background_counts - expected_model_counts - s2);
  return b;
  
  
}


vector pgstat(vector observed_counts, vector background_counts, vector background_error, vector expected_model_counts, int[] idx_background_zero, int[] idx_background_nonzero) {
  
  int N = num_elements(expected_model_counts);
  vector[N] log_likes;    
  vector[N] s2 = square(background_error);
  vector[N] b = background_model(observed_counts, background_counts, background_error, expected_model_counts);
  vector[N] factorial_term = expected_model_counts - lgamma(observed_counts +1 );
  
  log_likes[idx_background_nonzero] = (-square(b[idx_background_nonzero] - background_counts[idx_background_nonzero]) ./ (2 * s2[idx_background_nonzero])
				       + observed_counts[idx_background_nonzero] .* log(b[idx_background_nonzero] + expected_model_counts[idx_background_nonzero])
				       - b[idx_background_nonzero] - factorial_term[idx_background_nonzero]
				       - 0.5 * log(2 * pi()) - log(background_error[idx_background_nonzero]));	
  
  
  for (n in 1:num_elements(idx_background_zero)) {    
    log_likes[idx_background_zero[n]] = lmultiply(observed_counts[idx_background_zero[n]], expected_model_counts[idx_background_zero[n]]) - factorial_term[idx_background_zero[n]];	
  }
  
  return log_likes;
  
}
