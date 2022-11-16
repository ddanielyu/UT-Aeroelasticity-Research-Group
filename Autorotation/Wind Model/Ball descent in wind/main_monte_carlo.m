%This code performs Markov Chain Monte Carlo simulations of a ball falling from
%100 m in the presence of wind and gravity.
%The drag area is 0.14 m2 and initial velocity at 10 m height is 5 mph (2.2m/s).
%The hourly wind velocity data is obtained from Austin airport data at 10 m.
%A statistical profile at 10 Hz is obtained using Markov Chain method and
%extrapolated to other altitudes using log-law atmospheric wind profile.
%The equations of drag & gravity are simultaneously solved to obtain trajectory.
%Ref [1]: Papaefthymiou, G. and Klockl, B., 2008. MCMC for wind power simulation

%% Simulation Parameters
dropHeight = 100;%m
nSimulations = 100;%number of Monte Carlo simulations

%% Initializations
finalPosition = zeros(100,1);
trajectorySimulation = zeros(nSimulations,2,400);
lengthSimulation = zeros(nSimulations,1);

%% Read Wind Velocity Data
weatherData = readtable('weather_data.csv','Delimiter',',');
windData = zeros(size(weatherData,1),1);
for iWeatherData = 1:size(weatherData,1)
    windData(iWeatherData) = str2double(weatherData.Wind_Speed{iWeatherData}...
        (1:(end-4)));
end

%% Transition probability and cumulative probability matrix
% This section converts the wind data into state space form and generate
% the cumulative probability matrix as in Ref. [1]
transitionMatrix = [windData(1:end-1),windData(2:end)];
eventMatrix = zeros(max(windData)+1);
for iWeatherData = 1:length(windData)-1
    eventMatrix(transitionMatrix(iWeatherData,1)+1,...
        transitionMatrix(iWeatherData,2)+1) = ...
        eventMatrix(transitionMatrix(iWeatherData,1)+1,...
        transitionMatrix(iWeatherData,2)+1) + 1;
end
probabilityMatrix = eventMatrix.*repmat(1./sum(eventMatrix,2),1,...
    max(windData)+1);
probabilityMatrix = fillmissing(probabilityMatrix,'constant',0);
probabilityMatrixCumulative = cumsum(probabilityMatrix,2);


%% Trajectory vectors for each simulation and the final position at touchdown
for iSimulation = 1:nSimulations
    fprintf('Simulation %u of %u\n',iSimulation,nSimulations);
    trajectoryVector = model_dynamics(dropHeight,probabilityMatrixCumulative);
    finalPosition(iSimulation) = ...
        min(trajectoryVector(1,trajectoryVector(2,:)<0));
    lengthSimulation(iSimulation) = length(trajectoryVector);
    trajectorySimulation(iSimulation,:,1:lengthSimulation(iSimulation)) = ...
        real(trajectoryVector);
    
end

%% Plot of trajectories, mean & standard deviation. Histogram of final positions
figure(1)
hold all;
for iSimulation = 1:nSimulations
    plot(-squeeze(trajectorySimulation(iSimulation,1,...
        1:lengthSimulation(iSimulation))),...
        squeeze(trajectorySimulation(iSimulation,2,...
        1:lengthSimulation(iSimulation))),...
        '.','MarkerSize',1,'Color',[180,180,180]/255)
    
    meanTrajectorySimulation = ...
        mean(trajectorySimulation(:,:,1:lengthSimulation(iSimulation)),1);
    stdTrajectorySimulation = ...
        std(trajectorySimulation(:,:,1:lengthSimulation(iSimulation)),1);
    stdTrajectorySimulation(1,1,stdTrajectorySimulation(1,1,:)>15)=15;
    %removes outliers
end

trajectoryPlotHandle = plot(-squeeze(meanTrajectorySimulation(1,1,:)),...
    squeeze(meanTrajectorySimulation(1,2,:)),...
    '-k');

stdLeftHandle = plot(...
    -squeeze(meanTrajectorySimulation(1,1,:)-stdTrajectorySimulation(1,1,:)),...
    squeeze(meanTrajectorySimulation(1,2,:)),'--k');

stdRightHandle = plot(...
    -squeeze(meanTrajectorySimulation(1,1,:)+stdTrajectorySimulation(1,1,:)),...
    squeeze(meanTrajectorySimulation(1,2,:)),'--k');

legend([trajectoryPlotHandle,stdLeftHandle],'mean','1 standard deviation')
xlabel('Horizontal displacement [m]')
ylabel('Height above ground [m]')
title({['Monte Carlo simulation of ball dropped from ', ...
    num2str(dropHeight),' m'],['in 5 mph (=2.24 m/s) wind (A x C_D = 0.14)']})
axis equal
xlim([-10,80])
ylim([0,dropHeight])
box on

figure(2);
histogram(-finalPosition,20,'Normalization','probability')
xlabel('Horizontal displacement at touchdown [m]')
ylabel('Probability of Occurance [m]')
title({'Probability distribution of horizontal displacements in m',...
    ['(100 drops from ',num2str(dropHeight),' m in 5 mph (2.24 m/s) wind)']})