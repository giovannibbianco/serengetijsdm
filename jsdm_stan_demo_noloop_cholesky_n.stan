// following suggestions by nhuurre in Stan forum

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
  real<lower=0> area[n_sites];     // area for every site
  int<lower=1> K;                 // number of sample-level predictors
  int<lower=1> n_s;               // num of species
  int<lower=1> n_t;               // num species level predictors (traits)
  int<lower=1,upper=n_s> sp[n_obs];   // species id 
  int<lower=1,upper=n_sites> site[n_obs];
  matrix[n_sites, K] X;                 // obs-level design matrix 
  matrix[n_s, n_t] TT;            // species-level traits
  matrix[n_s, n_s] C;             // phylogenetic correlation matrix
  // vector[n_s] ones;               // vector on 1s
  int<lower=1> n_max[n_s];        // Upper bound of population size per spp
  real<lower=0,upper=1> p_obs[n_s];
}

transformed data {
  
  int<lower=0> Y[n_sites, n_s]; // total spp by site
  
  for(i in 1:n_sites){
    for(j in 1:n_s){
      Y[i,j] = 0;
    }
  }
  
  for(i in 1:n_obs){
    Y[site[i], sp[i]] += 1;
  }
}

parameters {
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0>[K] tau;

  matrix[n_t, K] Z;
  vector[n_s * K] betas;
  real<lower=0,upper=1> rho;      // correlation between phylogeny and betas
  real<lower=0,upper=1> p[n_s];   // detection probability
}

transformed parameters { 
  matrix[K, K] L_Sigma = diag_pre_multiply(tau,L_Omega);
  matrix[n_s*K, n_s*K] L_S = kronecker(L_Sigma,
          cholesky_decompose(rho * C + diag_matrix(rep_vector(1-rho, n_s))) );
          
  //matrix[K, K] L_Sigma = diag_pre_multiply(tau, L_Omega);
  //matrix[n_s, n_s] L_phylo = cholesky_decompose(
  //           rho * C + diag_matrix(rep_vector(1 - rho, K)));
 
 // matrix[n_t, K] Z = to_matrix(z, n_t, K);    
  //matrix[n_s,  K] m = TT * Z;        // mean of coeffs
  vector[n_s * K] m = to_vector(TT * Z); 
  matrix[n_s, K] b_m = to_matrix(betas, n_s, K);  // coeffs
  //matrix[n_s, K] beta = m + L_phylo * beta_std * L_Sigma'; // Matrix Normal NCP
} 

model {

  matrix[n_sites, n_s] log_lambda;
  int Ymax[n_sites, n_s];
  // priors
  // p ~ beta(2,2);
  L_Omega ~ lkj_corr_cholesky(2);
  tau ~ student_t(3,0,10); 
  betas ~ multi_normal_cholesky(m, L_S);
  //to_vector(beta_std) ~ std_normal();
  to_vector(Z) ~ normal(0, 2); 
  
  //rho ~ beta(2,2);
  
  // mix prior on rho
  //target += log_sum_exp(log(0.5) +  beta_lpdf(rho|1, 10), log(0.5) +  beta_lpdf(rho|2,2));
  
  for (n in 1:n_sites){
    for(s in 1:n_s){
      log_lambda[n,s] = dot_product( X[n,] , b_m[s,]) + log(area[n]);
      Ymax[n,s] = Y[n,s] + qpois(0.9999, exp(log_lambda[n,s]) * (1 - p_obs[s]), Y[n,s] + n_max[s]);
    }
  }
  
  for (n in 1:n_sites){
    for(s in 1:n_s){
      target += poisson_log_lpmf(Y[n,s] | log_lambda[n,s] + log(p_obs[s]));
    }
  }
}

generated quantities{
  int<lower=0> N[n_sites,n_s];
  real<lower=0> D[n_s];
  
  for (n in 1:n_sites){
    for(s in 1:n_s){
      N[n,s]= poisson_log_rng(dot_product( X[n,] , b_m[s,]) + log(area[n])); 
      }
    }
  
  for(s in 1:n_s) D[s] = sum(N[,s])/sum(area);
  
}