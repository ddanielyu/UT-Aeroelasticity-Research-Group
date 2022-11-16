function U = markov_chain(chain_len,init_vel)
%% Generate random numbers
load('wind_data.mat');
rand_nums = rand(chain_len,1);
%% Markov Chain generation
% init_vel = 6;%corresponds to wind speed of init_state-1 mph
U = (init_vel-1)*ones(chain_len+1,1);
for l = 2:chain_len
    l;
    U(l) =  find((rand_nums(l)-P_cumulative(U(l-1)+1,:))<=0,1)-1;
end
