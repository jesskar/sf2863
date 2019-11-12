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

%% Analytic Solution

%contains stationary distributions
Pi_matrix = []; 

for i=1:length(muvec_1)
    mu_1 = muvec_1(i);
    mu_2 = muvec_2(i);
    
    %choose ONE of below
    Q = bestQmatrix(mu_AB, mu_C);
    %Q = getQmatrix(mu_1,mu_2);
    
    P = expm(Q);
    [V,D] = eig(P');
    Pi = V(:,1)';
    Pi = Pi/sum(Pi);
    Pi_matrix(:,i) = Pi'; 
    
end

production_1 = d * Pi_matrix;

% Save data
format long g
Analytic = [Pi_matrix' production_1']
Analytic = round(Analytic,2)

%% Continuous time simulation inc. workers break

% converges with 1% for T = 100000
T_max = 100000; x = 0.6; 
vec = [1,2,3,4];

time_matrix = [];
probability = [];
for i=1:length(muvec_1)
    mu_1 = muvec_1(i);
    mu_2 = muvec_2(i);
   
    %choose ONE of below
    Q = bestQmatrix(mu_AB, mu_C);
    %Q = getQmatrix(mu_1,mu_2);
    
    t = 0; 
    state = datasample(vec,1);
    state_time = [0 0 0 0];
    both_working_visits = 0;
    longer_than_x = 0;
    
    while t < T_max
        
        T_min = Inf;
        next = 0;
        
        for tmp_state=1:length(Q)
            
            % We have to move. Look at neighbor states. Next.
            if tmp_state == state
                
                continue;
                
            end
                                    
            % One of the states has q_{ij}=0. We cannot move here. Next.
            if Q(state,tmp_state) == 0
                
                continue;
                
            end
            
            % Sample time
            T_i = exprnd(1/Q(state,tmp_state));
            
            if T_i < T_min
                
                T_min = T_i;
                next = tmp_state;
                
            end
            
            % Frequentist calc of prob longer than x
            if state == 1
            
                both_working_visits = both_working_visits + 1;

                if T_i > x

                    longer_than_x = longer_than_x + 1;

                end
            
            end
            
        end   
        
        % Update state
        state_time(state) = state_time(state) + T_min;
        t = t + T_min;
        state = next;
        
    end
    
    time_matrix(:,i) = (state_time/sum(state_time))';
    
    probability = [probability longer_than_x/both_working_visits];

    
end
   
production_2 = d * time_matrix;
%% Discretization-approach simulation

N=1000000;
h=0.001;
t_max = h*N;

state_matrix = [];
for i=1:length(muvec_1)
    mu_1 = muvec_1(i);
    mu_2 = muvec_2(i);
    
    %choose ONE of below
    %Q = bestQmatrix(mu_AB, mu_C);
    Q = getQmatrix(mu_1,mu_2);

    I = eye(4);
    P = I + Q*h; 
    
    vec = [1,2,3,4];
    state = datasample(vec,1);
    state_count = [0 0 0 0];
    t=0;
    
    while t < t_max
        
        state_count(state) = state_count(state) + 1;
        X = rand();
        
        prob_vec = [];
        for k=1:4
            prob_vec(k) = sum(P(state,vec(1:k)));
        end

        stats = find(prob_vec > X);
        state = stats(1);
        
        t=t+h;
        
    end
    
    state_matrix(:,i) = (state_count/sum(state_count))';

end

production_3 = d * state_matrix;