%We compute the intensity matrix Q (continuous time)
%Then translate it using P=exp(Qt) 

%clear all; clear workspace
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

%% Continuous time simulation
% converges with 1% for T = 100000
T_max = 100000;
time_matrix = [];
for i=1:length(muvec_1)
    mu_1 = muvec_1(i);
    mu_2 = muvec_2(i);
   
    %choose ONE of below
    Q = bestQmatrix(mu_AB, mu_C);
    %Q = getQmatrix(mu_1,mu_2);
    t = 0;
    state = 1; %Possible values, 1, 2, 3, 4.
    state_time = [0 0 0 0];
    
    while t < T_max
        
        T_min = Inf;
        next = 0;
        
        for tmp_state=1:length(Q)
            
            % We have to move. Look at neighbor states. Next.
            if tmp_state == state
                
                continue;
                
            end
                        
            Q(state,tmp_state)
            
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
            
        end   
        
        state_time(state) = state_time(state) + T_i;
        t = t + T_i;
        state = next;

    end
    % Sum over time in each state, not frequency
    time_matrix(:,i) = (state_time/sum(state_time))';
    
end
   
production_2A = d * time_matrix;