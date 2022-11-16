function velocityWindChainAt10 = generate_markov_chain(chainLength,...
    initialVelocity,probabilityMatrixCumulative)
%This function generates wind velocity at 10 m height at 10 Hz using Markov
%Chain method and measured wind velocity data at hourly intervals
%Ref [1]: Papaefthymiou, G. and Klockl, B., 2008. MCMC for wind power simulation
%% Generate random numbers
randomNumbers = rand(chainLength,1);
%% Markov Chain generation
velocityWindChainAt10 = (initialVelocity-1)*ones(chainLength+1,1);
for iChain = 2:chainLength
    velocityWindChainAt10(iChain) =  find((randomNumbers(iChain)-...
      probabilityMatrixCumulative(velocityWindChainAt10(iChain-1)+1,:))<=0,1)-1;
end
