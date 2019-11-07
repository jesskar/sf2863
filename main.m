%We compute the intensity matrix Q (continuous time)
%Then translate it using P=exp(Qt) t

% Parameters
rng(1)
d1=30; 
d2=90; 
d12=120; 
d_0=10; 
d = [d12 d2 d1 d_0];
T_max = 1000000; %h
state_init = 1; %Possible values, 1, 2, 3, 4.
h=0.001;

mu_A=0.70; 
mu_B=0.70; 
mu_C=0.50; 
mu_AB=1.75; 
mu_AC=1.30; 
mu_BC=1.3;

lambda_1=0.20; 
lambda_2=0.50; 
muvec_1 = [mu_AB mu_AC mu_C  mu_BC mu_A  mu_B mu_A];
muvec_2 = [mu_C  mu_B  mu_AB mu_A  mu_BC mu_AC mu_B];

Pimatrix = [];
prod2 = [];

for i=1:length(muvec_1)
    mu_1 = muvec_1(i);
    mu_2 = muvec_2(i);
    
    q00 = -lambda_1 - lambda_2 - (lambda_1 + lambda_2);
    q01 = lambda_1;
    q02 = lambda_2;
    q03 = (lambda_1 + lambda_2);

    q10 = mu_1;
    q11 = -mu_1 - lambda_2 - (mu_1 + lambda_2);
    q12 = (mu_1 + lambda_2);
    q13 = lambda_2;

    q20 = mu_2;
    q21 = (mu_2 + lambda_1);
    q22 = -mu_2 -lambda_1 - (mu_2 + lambda_1); 
    q23 = 1/lambda_1;

    q30 = (mu_1 + mu_2);
    q31 = mu_2;
    q32 = mu_1;
    q33 = -mu_1 - mu_2 - (mu_1 + mu_2);

    Q = [q00 q01 q02 q03;
         q10 q11 q12 q13;
         q20 q21 q22 q23;
         q30 q31 q32 q33];

    P = expm(Q);

    [V,D] = eig(P');
    Pi = V(:,1)';
    Pi = Pi./sum(Pi);
    Pimatrix(:,i) = Pi'; %Stationary distribution matrixprod
    
    state = state_init;
    state_count = [0 0 0 0];
    vec = [1,2,3,4];
    t = 0;
    tic
while t < T_max(end)
	vec(state) = [];
	X_1 = exprnd(Q(state,vec(1)));
	X_2 = exprnd(Q(state,vec(2)));
	X_3 = exprnd(Q(state,vec(3)));
	[X,state] = min([X_1 X_2 X_3]);
    
    state = vec(state);
    vec = [1,2,3,4];
    
    state_count(state) = state_count(state) + X;
    state_count = state_count/sum(state_count);
    prodcon = sum(d.*state_count);
	
    t = t + X;
end
toc
    prod2(i) = prodcon;
end

prod = d * Pimatrix;




