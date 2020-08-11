
functions { 
  
  /* compute the kronecker product
  * Args: 
  *   A,B: matrices 
  * Returns: 
  *   kronecker product of A and B
  */ 
  matrix kronecker(matrix A, matrix B) { 
    matrix[rows(A)*rows(B), cols(A)*cols(B)] kron; 
    for (i in 1:cols(A)) { 
      for (j in 1:rows(A)) { 
        kron[((j-1)*rows(B)+1):(j*rows(B)), ((i-1)*cols(B)+1):(i*cols(B))] = A[j,i] * B;
      } 
    } 
    return kron; 
  } 
  
  
  int qpois(real q, real lambda, int max_x) {
    int x = 0;
    real res = poisson_cdf(x, lambda);
    
    while(res < q && x < max_x){
      x = x + 1;
      res = poisson_cdf(x, lambda);
    }
    return x; 
  } 
} 

data { 
  int<lower=1> n_obs;
  int<lower=1> n_sites;           // total number of observations (sites/segments)
  int<lower=1> n_dates;
  int<lower=1> n_tcov;            // number of time varying covariates (2 for the moment)
  real<lower=0> area[n_sites];    // area for every site
  int<lower=1> K;                 // number of sample-level predictors
  int<lower=1> n_s;               // num of species
  int<lower=1> n_t;               // num species level predictors (traits)
  int<lower=1,upper=n_s> sp[n_obs];   // species id 
  int<lower=1,upper=n_sites> site[n_obs];
  matrix[n_sites, K] X;           // obs-level design matrix 
  matrix[n_sites, n_dates] Xt;
  matrix[n_s, n_t] TT;            // species-level traits
  matrix[n_s, n_s] C;             // phylogenetic correlation matrix
  vector[n_s] ones;               // vector on 1s
  int<lower=1> n_max[n_s];        // Upper bound of population size per spp
  real<lower=0,upper=1> p_obs[n_s];
  int<lower=1> date[n_obs];
  
}

transformed data {
  
  int<lower=0> Y[n_sites, n_s, n_dates]; // total spp by site and date
  
  for(i in 1:n_sites){
    for(j in 1:n_s){
      for(d in 1:n_dates){
        Y[i,j,d] = 0;
      }
    }
  }
  
  for(i in 1:n_obs){
    Y[site[i], sp[i], date[i]] += 1;
  }
}


parameters {
  corr_matrix[K+n_tcov] Omega;           // correlation matrix for var-covar of betas
  vector<lower=0>[K+n_tcov] tau;         // scales for the variance covariance of betas
  vector[n_s * (K+n_tcov)] betas;
  real<lower=0,upper=1> rho;             // correlation between phylogeny and betas
  vector[n_t * (K+n_tcov)] z;            // coeffs for traits
  real<lower=0,upper=1> p[n_s];          // detection probability
}

transformed parameters { 
  matrix[K+n_tcov, K+n_tcov] Sigma = quad_form_diag(Omega, tau);
  matrix[n_s*(K+n_tcov),n_s*(K+n_tcov)] S=kronecker(Sigma, rho*C+(1-rho)*diag_matrix(ones));
  matrix[n_t, K+n_tcov] Z = to_matrix(z, n_t, K+n_tcov);    
  vector[n_s * (K+n_tcov)] m = to_vector(TT * Z);        // mean of coeffs
  matrix[n_s, K+n_tcov] b_m = to_matrix(betas, n_s, K+n_tcov);  // coeffs
} 

model {
  real log_lambda[n_sites, n_s, n_dates];
  int Ymax[n_sites, n_s, n_dates];
  // priors
  // p ~ beta(2,2);
  Omega ~ lkj_corr(2);
  tau ~ student_t(3,0,10); // cauchy(0, 2.5); // lognormal()
  betas ~ multi_normal(m, S);
  //rho ~ beta(2,2);
  z ~ normal(0,2);
  
  // mix prior on rho
  //target += log_sum_exp(log(0.5) +  beta_lpdf(rho|1, 10), log(0.5) +  beta_lpdf(rho|2,2));
  
  for (n in 1:n_sites){
    for(s in 1:n_s){
      for(d in 1:n_dates){
        log_lambda[n,s,d] = dot_product( X[n,] , b_m[s,1:K]) 
        + Xt[n,d]*b_m[s,K+n_tcov] 
        + log(area[n]);
        Ymax[n,s,d] = Y[n,s,d] + qpois(0.999, exp(log_lambda[n,s,d])*(1-p_obs[s]), Y[n,s,d]+n_max[s]);
      }
    }
  }
  
  for (n in 1:n_sites){
    for(s in 1:n_s){
      for(d in 1:n_dates){
        vector[Ymax[n,s,d] - Y[n,s,d] + 1] lp;
        for (j in 1:(Ymax[n,s,d]  - Y[n,s,d] + 1)){
          lp[j] = poisson_log_lpmf(Y[n,s,d] + j - 1 | log_lambda[n,s,d]) 
          + binomial_lpmf(Y[n,s,d] | Y[n,s,d] + j - 1, p_obs[s]);
          }
          target += log_sum_exp(lp);
      }
    }
  }
}

generated quantities{
  int<lower=0> N[n_sites, n_s, n_dates]; 
  real<lower=0> D[n_s, n_dates];

  for (n in 1:n_sites){
    for(s in 1:n_s){
      for(d in 1:n_dates){
        N[n,s,d]= poisson_log_rng(dot_product( X[n,] , b_m[s,1:K]) 
        + Xt[n,d]*b_m[s,K+n_tcov] 
        + log(area[n]));
      }
    }

  for(s in 1:n_s){
    for(d in 1:n_dates){
      D[s,d] = sum(N[,s,d])/sum(area);
    }
  } 

}

