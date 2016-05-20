
/*
## ======================================================================
## Copyright 2004-2016, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================
*/


/*
## ======================================================================
From http://cran.r-project.org/doc/manuals/r-release/R-exts.html:

"Note that the exponential and gamma distributions are parametrized by
 scale rather than rate."
======================================================================
 */

#include <R.h>
#include <Rmath.h>


double threep_pmingamma_over_p (double x, int exposures, double *p,  double alpha,  double beta, 
				int lower_tail, int give_log) {
  // the cdf of the minimum gamma marginalized over p.
  // if lower_tail is non-zero calculate one minus the cdf instead.
  // if give_log is non-zero calculate the log density instead.
  // See the file 'Minimum_gamma_pdf_and_cdf.pdf' for further details.
  
  double p1, p2, p3, q1, q2;
  double G, G2, G3, G4;
  double one_minus_G, one_minus_G2, one_minus_G3, one_minus_G4;

  p1 = p[0];
  p2 = p[1];
  p3 = p[2];
  q1 = 1.0 - p1;
  q2 = 1.0 - p2;

  one_minus_G  = pgamma(x, alpha, beta, FALSE, FALSE);
  one_minus_G2 = one_minus_G  * one_minus_G;
  one_minus_G3 = one_minus_G2 * one_minus_G;
  one_minus_G4 = one_minus_G3 * one_minus_G;

  G  = 1.0 - one_minus_G;
  G2 = 1.0 - one_minus_G2;
  G3 = 1.0 - one_minus_G3;
  G4 = 1.0 - one_minus_G4;

  if (exposures == 2) {

    G = G * q1 + G2 * p1;

  } else if (exposures == 3) {

    G = G * q1 * q1 + G2 * p1 * (q1 + q2) + G3 * p1 * p2;

  } else if (exposures == 4) {

    G = G * q1 * q1 * q1 +
      G2 * p1 * ( q2 * q2 + q1 * q2 + q1 * q1 ) +
      G3 * p1 * p2 * ( q1 + q2 + 1.0 - p3 ) +
      G4 * p1 * p2 * p3;
  }

  if (lower_tail) {

    if (give_log)
      return log(G);
    else
      return G;

  } else {

    if (give_log) 
      return log(1.0-G);
    else
      return (1.0-G);

  }
}



double threep_dmingamma_over_p (double x, int exposures, double *p,  
				double alpha,  double beta, int give_log) {

  // density of the minimum gamma marginalized over p.
  // if give_log is non-zero calculate the log density instead.
  // See the file 'Minimum_gamma_pdf_and_cdf.pdf' for further details.

  double dens;
  double g;
  double p1, p2, p3, q1, q2;
  double g2, g3, g4, one_minus_F;

  p1 = p[0];
  p2 = p[1];
  p3 = p[2];
  q1 = 1.0 - p1;
  q2 = 1.0 - p2;

  g = dgamma(x, alpha, beta, FALSE);
  one_minus_F = pgamma(x, alpha, beta, FALSE, FALSE);

  g2 = g * 2.0 * one_minus_F;
  g3 = g * 3.0 * one_minus_F * one_minus_F;
  g4 = g * 4.0 * one_minus_F * one_minus_F * one_minus_F;

  if (exposures == 1) {

    dens = g;

  } else if (exposures == 2) {
    
    dens = g * q1 + g2 * p1;

  } else if (exposures == 3) {

    dens = g * q1 * q1 + g2 * p1 * (q1 + q2) + g3 * p1 * p2;

  } else if (exposures == 4) {

    dens = g * q1 * q1 * q1 +
      g2 * p1 * ( q2 * q2 + q1 * q2 + q1 * q1 ) +
      g3 * p1 * p2 * ( q1 + q2 + 1.0 - p3 ) +
      g4 * p1 * p2 * p3;
  }

  if (give_log) 
    return log(dens);
  else
    return dens;
}




void threep_RT_log_likelihood (double *loglik, int *n,
			       double *RT, double *T0,
			       int *resp_new, int *exposures,
			       double *p, double *alpha, double *beta) {
  
  // See the RT_likelihood.R file for further details of the arguments.
  
  // beta = ( Oold, Nold, Onew, Nnew )
  
  int t;
  double z, alpha_t, exposures_t;
  
  for (t = 0; t < (*n); t++) {
    
    z = RT[t] - (*T0);
    alpha_t = alpha[t];
    exposures_t = exposures[t];
    
    if (exposures_t > 0) { // old stimulus?
      
      if (resp_new[t] == 0) { // respond old
	
	loglik[t] = threep_dmingamma_over_p(z, exposures_t, p, alpha_t, beta[0], TRUE) + 
	  pgamma(z, alpha_t, beta[1], FALSE, TRUE);
      	
      } else { // respond new
	
	loglik[t] = dgamma(z, alpha_t, beta[1], TRUE) + 
	  threep_pmingamma_over_p(z, exposures_t, p, alpha_t, beta[0], FALSE, TRUE);
      }
      
    } else { // new stimulus?
      
      if (resp_new[t] == 0) { // respond old
	
	loglik[t] = dgamma(z, alpha_t, beta[2], TRUE) + 
	  pgamma(z, alpha_t, beta[3], FALSE, TRUE);
	
      } else { // respond new
	
	loglik[t] = dgamma(z, alpha_t, beta[3], TRUE) + 
	  pgamma(z, alpha_t, beta[2], FALSE, TRUE);
	
      }
    }
  }
}



double fourp_pmingamma_over_p (double x, int exposures, double *p,  double alpha,  double beta, double beta_star,
			       int lower_tail, int give_log) {
  
  // the cdf of the minimum gamma marginalized over p.
  // if lower_tail is non-zero calculate one minus the cdf instead.
  // if give_log is non-zero calculate the log density instead.
  // See the file 'Minimum_gamma_pdf_and_cdf.pdf' for further details.
  
  double p0, p1, p2, p3, q0, q1, q2, q3;
  double G0, G1, G2, G3, G4;
  double one_minus_G, one_minus_G_star;
  double the_cdf;

  p0 = p[0];
  p1 = p[1];
  p2 = p[2];
  p3 = p[3];
  q0 = 1.0 - p0;
  q1 = 1.0 - p1;
  q2 = 1.0 - p2;
  q3 = 1.0 - p3;
  
  one_minus_G_star  = pgamma(x, alpha, beta_star, FALSE, FALSE);  
  one_minus_G       = pgamma(x, alpha, beta, FALSE, FALSE);

  G0 = 1.0 - one_minus_G_star;
  G1 = 1.0 - one_minus_G_star * one_minus_G;
  G2 = 1.0 - one_minus_G_star * one_minus_G * one_minus_G;
  G3 = 1.0 - one_minus_G_star * one_minus_G * one_minus_G * one_minus_G;
  G4 = 1.0 - one_minus_G_star * one_minus_G * one_minus_G * one_minus_G * one_minus_G;
  
  if (exposures == 1) {

    the_cdf = G0 * q0 +
      G1 * p0;
    
  } else if (exposures == 2) {

    the_cdf = G0 * q0 * q0 +
      G1 * p0 * (q0 + q1) +
      G2 * p0 * p1;

  } else if (exposures == 3) {

    the_cdf = G0 * q0 * q0 * q0 +
      G1 * p0 * ( q1 * q1 + q0 * q1 + q0 * q0 ) +
      G2 * p0 * p1 * ( q0 + q1 + q2 ) +
      G3 * p0 * p1 * p2;

  } else if (exposures == 4) {

    the_cdf =  G0 * (q0 * q0 * q0 * q0) +
      G1 * p0 * ( (q0 * q0 * q0) + (q0 * q0 * q1) + (q0 * q1 * q1) + (q1 * q1 * q1) ) +
      G2 * p0 * p1 * ( (q0 * q0) + (q0 * q1) + (q0 * q2) + (q1 * q1) + (q1 * q2) + (q2 * q2) ) +
      G3 * p0 * p1 * p2 * (q0 + q1 + q2 + q3) +
      G4 * p0 * p1 * p2 * p3;  
  }

  if (lower_tail) {

    if (give_log)
      return log(the_cdf);
    else
      return the_cdf;

  } else {

    if (give_log) 
      return log(1.0-the_cdf);
    else
      return (1.0-the_cdf);

  }
}



double fourp_dmingamma_over_p (double x, int exposures, double *p,  
			       double alpha,  double beta, double beta_star,
			       int give_log) {

  // density of the minimum gamma marginalized over p.
  // if give_log is non-zero calculate the log density instead.
  // See the file 'Minimum_gamma_pdf_and_cdf.pdf' for further details.

  double the_pdf;
  double p0, p1, p2, p3, q0, q1, q2, q3;
  double g, g_star;
  double g0, g1, g2, g3, g4;
  double one_minus_G_star, one_minus_G, one_minus_G2, one_minus_G3, one_minus_G4;

  p0 = p[0];
  p1 = p[1];
  p2 = p[2];
  p3 = p[3];
  q0 = 1.0 - p0;
  q1 = 1.0 - p1;
  q2 = 1.0 - p2;
  q3 = 1.0 - p3;

  g      = dgamma(x, alpha, beta,      FALSE);
  g_star = dgamma(x, alpha, beta_star, FALSE);
  
  one_minus_G_star  = pgamma(x, alpha, beta_star, FALSE, FALSE);  
  one_minus_G       = pgamma(x, alpha, beta, FALSE, FALSE);

  one_minus_G2 = one_minus_G  * one_minus_G;
  one_minus_G3 = one_minus_G2 * one_minus_G;
  one_minus_G4 = one_minus_G3 * one_minus_G;

  g0 = g_star;
  g1 = g_star * one_minus_G  + 1.0 * g * one_minus_G_star;
  g2 = g_star * one_minus_G2 + 2.0 * g * one_minus_G_star * one_minus_G;
  g3 = g_star * one_minus_G3 + 3.0 * g * one_minus_G_star * one_minus_G2;
  g4 = g_star * one_minus_G4 + 4.0 * g * one_minus_G_star * one_minus_G3;

  if (exposures == 1) {

    the_pdf = g0 * q0 +
      g1 * p0;
    
  } else if (exposures == 2) {

    the_pdf = g0 * q0 * q0 +
      g1 * p0 * (q0 + q1) +
      g2 * p0 * p1;

  } else if (exposures == 3) {

    the_pdf = g0 * q0 * q0 * q0 +
      g1 * p0 * ( q1 * q1 + q0 * q1 + q0 * q0 ) +
      g2 * p0 * p1 * ( q0 + q1 + q2 ) +
      g3 * p0 * p1 * p2;

  } else if (exposures == 4) {

    the_pdf =  g0 * (q0 * q0 * q0 * q0) +
      g1 * p0 * ( (q0 * q0 * q0) + (q0 * q0 * q1) + (q0 * q1 * q1) + (q1 * q1 * q1) ) +
      g2 * p0 * p1 * ( (q0 * q0) + (q0 * q1) + (q0 * q2) + (q1 * q1) + (q1 * q2) + (q2 * q2) ) +
      g3 * p0 * p1 * p2 * (q0 + q1 + q2 + q3) +
      g4 * p0 * p1 * p2 * p3;  
  }

  if (give_log) 
    return log(the_pdf);
  else
    return the_pdf;
}




void fourp_RT_log_likelihood (double *loglik, int *n,
			      double *RT, double *T0,
			      int *resp_new, int *exposures,
			      double *p, double *alpha, double *beta) {
  
  // See the RT_likelihood.R file for further details of the arguments.

  // beta = ( Oold, Nold, Onew, Nnew, beta_star )

  int t;
  double z, alpha_t, exposures_t;

  for (t = 0; t < (*n); t++) {

    z = RT[t] - (*T0);
    alpha_t = alpha[t];
    exposures_t = exposures[t];

    if (exposures_t > 0) { // old stimulus?

      if (resp_new[t] == 0) { // respond old
	
	loglik[t] = fourp_dmingamma_over_p(z, exposures_t, p, alpha_t, beta[0], beta[4], TRUE) + 
	  pgamma(z, alpha_t, beta[1], FALSE, TRUE);
      	
      } else { // respond new
	
	loglik[t] = dgamma(z, alpha_t, beta[1], TRUE) + 
	  fourp_pmingamma_over_p(z, exposures_t, p, alpha_t, beta[0], beta[4], FALSE, TRUE);
      }
      
    } else { // new stimulus?
      
      if (resp_new[t] == 0) { // respond old
	
	loglik[t] = dgamma(z, alpha_t, beta[2], TRUE) + 
	  pgamma(z, alpha_t, beta[3], FALSE, TRUE);
	
      } else { // respond new
	
	loglik[t] = dgamma(z, alpha_t, beta[3], TRUE) + 
	  pgamma(z, alpha_t, beta[2], FALSE, TRUE);
	
      }
    }
  }
}




double fourp_pminlnorm_over_p (double x, int exposures, double *p,  double log_alpha,  double beta, double beta_star,
			       int lower_tail, int give_log) {
  
  // the cdf of the minimum lnorm marginalized over p.
  // if lower_tail is non-zero calculate one minus the cdf instead.
  // if give_log is non-zero calculate the log density instead.
  // See the file 'Minimum_lnorm_pdf_and_cdf.pdf' for further details.
  
  double p0, p1, p2, p3, q0, q1, q2, q3;
  double G0, G1, G2, G3, G4;
  double one_minus_G, one_minus_G_star;
  double the_cdf;

  p0 = p[0];
  p1 = p[1];
  p2 = p[2];
  p3 = p[3];
  q0 = 1.0 - p0;
  q1 = 1.0 - p1;
  q2 = 1.0 - p2;
  q3 = 1.0 - p3;
  
  one_minus_G_star  = plnorm(x, log_alpha, beta_star, FALSE, FALSE);  
  one_minus_G       = plnorm(x, log_alpha, beta, FALSE, FALSE);

  G0 = 1.0 - one_minus_G_star;
  G1 = 1.0 - one_minus_G_star * one_minus_G;
  G2 = 1.0 - one_minus_G_star * one_minus_G * one_minus_G;
  G3 = 1.0 - one_minus_G_star * one_minus_G * one_minus_G * one_minus_G;
  G4 = 1.0 - one_minus_G_star * one_minus_G * one_minus_G * one_minus_G * one_minus_G;
  
  if (exposures == 1) {

    the_cdf = G0 * q0 +
      G1 * p0;
    
  } else if (exposures == 2) {

    the_cdf = G0 * q0 * q0 +
      G1 * p0 * (q0 + q1) +
      G2 * p0 * p1;

  } else if (exposures == 3) {

    the_cdf = G0 * q0 * q0 * q0 +
      G1 * p0 * ( q1 * q1 + q0 * q1 + q0 * q0 ) +
      G2 * p0 * p1 * ( q0 + q1 + q2 ) +
      G3 * p0 * p1 * p2;

  } else if (exposures == 4) {

    the_cdf =  G0 * (q0 * q0 * q0 * q0) +
      G1 * p0 * ( (q0 * q0 * q0) + (q0 * q0 * q1) + (q0 * q1 * q1) + (q1 * q1 * q1) ) +
      G2 * p0 * p1 * ( (q0 * q0) + (q0 * q1) + (q0 * q2) + (q1 * q1) + (q1 * q2) + (q2 * q2) ) +
      G3 * p0 * p1 * p2 * (q0 + q1 + q2 + q3) +
      G4 * p0 * p1 * p2 * p3;  
  }

  if (lower_tail) {

    if (give_log)
      return log(the_cdf);
    else
      return the_cdf;

  } else {

    if (give_log) 
      return log(1.0-the_cdf);
    else
      return (1.0-the_cdf);

  }
}



double fourp_dminlnorm_over_p (double x, int exposures, double *p,  
			       double log_alpha,  double beta, double beta_star,
			       int give_log) {

  // density of the minimum lnorm marginalized over p.
  // if give_log is non-zero calculate the log density instead.
  // See the file 'Minimum_lnorm_pdf_and_cdf.pdf' for further details.

  double the_pdf;
  double p0, p1, p2, p3, q0, q1, q2, q3;
  double g, g_star;
  double g0, g1, g2, g3, g4;
  double one_minus_G_star, one_minus_G, one_minus_G2, one_minus_G3, one_minus_G4;

  p0 = p[0];
  p1 = p[1];
  p2 = p[2];
  p3 = p[3];
  q0 = 1.0 - p0;
  q1 = 1.0 - p1;
  q2 = 1.0 - p2;
  q3 = 1.0 - p3;

  g      = dlnorm(x, log_alpha, beta,      FALSE);
  g_star = dlnorm(x, log_alpha, beta_star, FALSE);
  
  one_minus_G_star  = plnorm(x, log_alpha, beta_star, FALSE, FALSE);  
  one_minus_G       = plnorm(x, log_alpha, beta, FALSE, FALSE);

  one_minus_G2 = one_minus_G  * one_minus_G;
  one_minus_G3 = one_minus_G2 * one_minus_G;
  one_minus_G4 = one_minus_G3 * one_minus_G;

  g0 = g_star;
  g1 = g_star * one_minus_G  + 1.0 * g * one_minus_G_star;
  g2 = g_star * one_minus_G2 + 2.0 * g * one_minus_G_star * one_minus_G;
  g3 = g_star * one_minus_G3 + 3.0 * g * one_minus_G_star * one_minus_G2;
  g4 = g_star * one_minus_G4 + 4.0 * g * one_minus_G_star * one_minus_G3;

  if (exposures == 1) {

    the_pdf = g0 * q0 +
      g1 * p0;
    
  } else if (exposures == 2) {

    the_pdf = g0 * q0 * q0 +
      g1 * p0 * (q0 + q1) +
      g2 * p0 * p1;

  } else if (exposures == 3) {

    the_pdf = g0 * q0 * q0 * q0 +
      g1 * p0 * ( q1 * q1 + q0 * q1 + q0 * q0 ) +
      g2 * p0 * p1 * ( q0 + q1 + q2 ) +
      g3 * p0 * p1 * p2;

  } else if (exposures == 4) {

    the_pdf =  g0 * (q0 * q0 * q0 * q0) +
      g1 * p0 * ( (q0 * q0 * q0) + (q0 * q0 * q1) + (q0 * q1 * q1) + (q1 * q1 * q1) ) +
      g2 * p0 * p1 * ( (q0 * q0) + (q0 * q1) + (q0 * q2) + (q1 * q1) + (q1 * q2) + (q2 * q2) ) +
      g3 * p0 * p1 * p2 * (q0 + q1 + q2 + q3) +
      g4 * p0 * p1 * p2 * p3;  
  }

  if (give_log) 
    return log(the_pdf);
  else
    return the_pdf;
}


void fourp_RT_log_likelihood_L (double *loglik, int *n,
				double *RT, double *T0,
				int *resp_new, int *exposures,
				double *p, double *alpha, double *beta) {
  
  // See the RT_likelihood.R file for further details of the arguments.

  // beta = ( Oold, Nold, Onew, Nnew, beta_star )

  int t;
  double z, log_alpha_t, exposures_t;

  for (t = 0; t < (*n); t++) {

    z = RT[t] - (*T0);
    log_alpha_t = log(alpha[t]);
    exposures_t = exposures[t];

    if (exposures_t > 0) { // old stimulus?

      if (resp_new[t] == 0) { // respond old
	
	loglik[t] = fourp_dminlnorm_over_p(z, exposures_t, p, log_alpha_t, sqrt(beta[0]), sqrt(beta[4]), TRUE) + 
	  plnorm(z, log_alpha_t, sqrt(beta[1]), FALSE, TRUE);
      	
      } else { // respond new
	
	loglik[t] = dlnorm(z, log_alpha_t, sqrt(beta[1]), TRUE) + 
	  fourp_pminlnorm_over_p(z, exposures_t, p, log_alpha_t, sqrt(beta[0]), sqrt(beta[4]), FALSE, TRUE);
      }
      
    } else { // new stimulus?
      
      if (resp_new[t] == 0) { // respond old
	
	loglik[t] = dlnorm(z, log_alpha_t, sqrt(beta[2]), TRUE) + 
	  plnorm(z, log_alpha_t, sqrt(beta[3]), FALSE, TRUE);
	
      } else { // respond new
	
	loglik[t] = dlnorm(z, log_alpha_t, sqrt(beta[3]), TRUE) + 
	  plnorm(z, log_alpha_t, sqrt(beta[2]), FALSE, TRUE);
	
      }
    }
  }
}


		


void Weibull_RT_log_likelihood (double *loglik, int *n,
				double *RT, double *T0,
				int *resp_new, int *exposures,
				double *alpha, double *beta) {
    
  int t;
  double z, alpha_t, exposures_t, old_scale;
  
  for (t = 0; t < (*n); t++) {
    
    z = RT[t] - (*T0);
    alpha_t = alpha[t];
    exposures_t = exposures[t];

    old_scale = beta[1] * pow(exposures_t, -1.0 * beta[0]);
    
    if (exposures_t > 0) { // old stimulus?
      
      if (resp_new[t] == 0) { // respond old
	
	loglik[t] = dweibull(z, alpha_t, old_scale, TRUE) +
	  pweibull(z, alpha_t, beta[2], FALSE, TRUE);
      	
      } else { // respond new
	
	loglik[t] = dweibull(z, alpha_t, beta[2], TRUE) + 
	  pweibull(z, alpha_t, old_scale, FALSE, TRUE);
      }
      
    } else { // new stimulus?
      
      if (resp_new[t] == 0) { // respond old
	
	loglik[t] = dweibull(z, alpha_t, beta[3], TRUE) + 
	  pweibull(z, alpha_t, beta[4], FALSE, TRUE);
	
      } else { // respond new
	
	loglik[t] = dweibull(z, alpha_t, beta[4], TRUE) + 
	  pweibull(z, alpha_t, beta[3], FALSE, TRUE);
	
      }
    }
  }
}





void count_xis (int *xi,
		int *num_trials,
		int *block_sel,
		int *num_blocks,
		int *xi0,
		int *xi1,
		int *xi2) {

  int b, t;
  
  for (b=0; b<(*num_blocks); b++) {

    xi0[b] = xi1[b] = xi2[b] = 0;

    for (t=0; t<(*num_trials); t++) {

      if (block_sel[t] == (b+1)) {

	if (xi[t]==0) xi0[b]++;
	if (xi[t]==1) xi1[b]++;
	if (xi[t]==2) xi2[b]++;
	
      }

    }
  }
}
		




void sample_q (int *xi,
	       int *num_trials,
	       int *block_sel,
	       int *num_blocks,
	       double *qprior,
	       double *q0,
	       double *q1,
	       double *q2) {

  int b, t;
  int xi0, xi1, xi2;
  double a0, a1, a2, sum_a;
  
  GetRNGstate();

  for (b=0; b<(*num_blocks); b++) {

    xi0 = xi1 = xi2 = 0;

    for (t=0; t<(*num_trials); t++) {

      if (block_sel[t] == (b+1)) {

	if (xi[t]==0) xi0++;
	if (xi[t]==1) xi1++;
	if (xi[t]==2) xi2++;

      }      
    }

    a0 = rgamma(qprior[0] + (double)xi0, 1.0);
    a1 = rgamma(qprior[1] + (double)xi1, 1.0);
    a2 = rgamma(qprior[2] + (double)xi2, 1.0);

    sum_a = a0 + a1 + a2;

    q0[b] = a0/sum_a;
    q1[b] = a1/sum_a;
    q2[b] = a2/sum_a;
  }

  PutRNGstate();  

}
		
