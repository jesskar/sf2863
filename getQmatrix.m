function Q = getQmatrix(mu_1, mu_2)
    global lambda_1 lambda_2
    q00 = -lambda_1 - lambda_2;
    q01 = lambda_1;
    q02 = lambda_2;
    q03 = 0;

    q10 = mu_1;
    q11 = -mu_1 - lambda_2;
    q12 = 0;
    q13 = lambda_2;

    q20 = mu_2;
    q21 = 0;
    q22 = -mu_2 -lambda_1; 
    q23 = lambda_1;

    q30 = 0;
    q31 = mu_2;
    q32 = mu_1;
    q33 = -mu_1 - mu_2;

    Q = [q00 q01 q02 q03;
         q10 q11 q12 q13;
         q20 q21 q22 q23;
         q30 q31 q32 q33];
end