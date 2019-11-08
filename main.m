%We compute the intensity matrix Q (continuous time)
%Then translate it using P=exp(Qt) 

clear all; clear workspace
% Parameters
global lambda_1 lambda_2
rng(1)
d1=30; 
d2=90; 
d12=120; 
d_0=10; 
d = [d12 d2 d1 d_0];
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

%%
Pi_matrix = []; %contains stationary distributions

for i=1:length(muvec_1)
    mu_1 = muvec_1(i);
    mu_2 = muvec_2(i);
    
    %choose ONE of below
    Q = bestQmatrix(mu_AB, mu_C);
    %Q = getQmatrix(mu_1,mu_2);
    
    P = expm(Q);

    [V,D] = eig(P');
    Pi = V(:,1)';
    Pi = Pi./sum(Pi);
    Pi_matrix(:,i) = Pi'; 
end

production_1 = d * Pi_matrix;

%% Continuous time simulation
% converges with 1% for T = 10000 vs. T = 50000
T_max = 50000;
time_matrix = [];
for i=1:length(muvec_1)
    mu_1 = muvec_1(i);
    mu_2 = muvec_2(i);
    state = 1; %Possible values, 1, 2, 3, 4.
    state_time = [0 0 0 0];
    vec = [1,2,3,4];
    
    %choose ONE of below
    Q = bestQmatrix(mu_AB, mu_C);
    Q = getQmatrix(mu_1,mu_2);

    t = 0;
    while t < T_max
        vec = [1,2,3,4];
        zero = find(Q(state,:)==0);
        vec([state zero]) = [];
        
        X_1 = exprnd(1/Q(state,vec(1)));
        X_2 = exprnd(1/Q(state,vec(2)));
        [X,tmpstate] = min([X_1 X_2]);
        
        %update state
        state = vec(tmpstate);
        
        %increment time
        state_time(state) = state_time(state) + X;
        t = t + 1;
    end
    
    % Sum over time in each state, not frequency
    time_matrix(:,i) = (state_time/sum(state_time))';
    
end
   
production_2 = d * time_matrix;

%% Discretization-approach simulation
N=100000;
h=0.001;
t_max = h*N;
prod2 = [];

for i=1:length(muvec_1)
    mu_1 = muvec_1(i);
    mu_2 = muvec_2(i);
    
    Q = getQmatrix(mu_1,mu_2);

    I = eye(4);
    P = I + Q*h;

    %sample initial state
    vec = [1,2,3,4];
    state = datasample(vec,1);
    state_count = [0 0 0 0];
    tic
    t=0;
    
    while t < t_max
        X = rand();
        %compare random sample to probabilities p_{ij}
        state = vec(state);
        prob_vec = [];
        for i=1:4
            prob_vec(i) = sum(P(state,vec(1:i)));
        end

        stats = find(prob_vec > X);
        i = stats(1);

        state_count(i) = state_count(i) + 1;
        %state_count = state_count/sum(state_count);
        prodcon = sum(d.*state_count);

        t=t+h;
    end

    prod2(i) = prodcon;
end

prod = d * P;

