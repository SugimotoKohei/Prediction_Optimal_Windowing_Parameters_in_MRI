data {
    int N;                // num of sample size
    int O;                // num of MRI operators
    int C;                // num of Imaging conditions
    int Idx_operator[N];  // indices of MRI operators
    int Idx_condition[N]; // indices of imaging conditions
    matrix[O,O] X;        // design matrix
    vector[2] Y[N];       // standardized windowing parameters (WL and WW)
}

parameters {
    vector[3] beta_WW;
    vector[3] beta_WL;
    matrix[C,O] r_WW;
    matrix[C,O] r_WL;
    vector<lower=0>[3] sigma_r_WW;
    vector<lower=0>[3] sigma_r_WL;
    vector<lower=0>[2] sigma_vec[C];
    cholesky_factor_corr[2] corr_chol[C];
}

transformed parameters {
    matrix[C,O] mu_WW;
    matrix[C,O] mu_WL;
    vector[2] mu[C,O];
    cholesky_factor_cov[2] cov_chol[C];

    for (c in 1:C) {
        cov_chol[c] = diag_pre_multiply(sigma_vec[c], corr_chol[c]);
        for (o in 1:O) {
            mu_WW[c,o] = dot_product(X[o], beta_WW) + r_WW[c,o];
            mu_WL[c,o] = dot_product(X[o], beta_WL) + r_WL[c,o];
            mu[c,o] = [mu_WW[c,o], mu_WL[c,o]]';
        }
    }
}

model {
    // Prior distribution
    for (c in 1:C) {
        target += lkj_corr_cholesky_lpdf(corr_chol[c] | 1);
        for (o in 1:O) {
            target += normal_lpdf(r_WW[c,o] | 0, sigma_r_WW[o]);
            target += normal_lpdf(r_WL[c,o] | 0, sigma_r_WL[o]);
        }
    }

    // Likelihood
    for (n in 1:N)
        target += multi_normal_cholesky_lpdf(Y[n] | mu[Idx_condition[n], Idx_operator[n]], cov_chol[Idx_condition[n]]);
}

generated quantities {
    matrix[2,2] corr[C];
    matrix[2,2] cov[C];
    vector[2] Y_pred[N];
    vector[2] Y_pred_operator[O,C];
    
    for (c in 1:C) {
        corr[c] = multiply_lower_tri_self_transpose(corr_chol[c]);
        cov[c] = multiply_lower_tri_self_transpose(cov_chol[c]);
    }

    for (n in 1:N) {
        Y_pred[n] = multi_normal_cholesky_rng(mu[Idx_condition[n], Idx_operator[n]], cov_chol[Idx_condition[n]]);
    }
    
    for (c in 1:C) {
        for (o in 1:O) {
            Y_pred_operator[o,c] = multi_normal_cholesky_rng(mu[c,o], cov_chol[c]);
        }
    }
}
