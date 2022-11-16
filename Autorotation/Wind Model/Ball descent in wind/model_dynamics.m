function positionBallVector = model_dynamics(dropHeight,...
    probabilityMatrixCumulative)
%This function solves the equations of horizontal and vertical forces on
%the ball to obtain the position vector as a function of time.

%% Constants, Simulation Parameters and Initializations
GRAVITY = 9.8;
DENSITY_AIR = 1.225;
heightReference = 10;
heightObstacles = 0.4;%https://websites.pmc.ucsc.edu/~jnoble/wind/extrap/
massBall = convmass(1.5, 'lbm','kg');
dragArea = 0.14;
velocityWindInitialAt10 = 5;%mph
dt = 0.1;
timeVector = 0:0.1:30;
nTimeSteps = length(timeVector);%length of markov chain: 30 s @ 10 Hz
velocityBallVector = zeros(2,nTimeSteps);
positionBallVector = repmat([0;dropHeight],1,nTimeSteps);%m
velocityBallInitial = [0;0];%m/s
jacobianMatrix = zeros(2);

%% Markov Chain wind velocities in m/s at 10 m
velocityWindChainAt10 = generate_markov_chain(nTimeSteps,...
    velocityWindInitialAt10+1,probabilityMatrixCumulative)*0.44704;

%% Time Marching Simulation
for iTimeStep = 2:nTimeSteps
    nIteration = 0;
    errorMagnitude = 1;
    while errorMagnitude > 1e-12 && nIteration < 30
        nIteration = nIteration+1;
        
        velocityWindAtH = velocityWindChainAt10(iTimeStep)*...
            log(positionBallVector(2,iTimeStep-1)/heightObstacles)...
            /log(heightReference/heightObstacles);
        %log-law wind atmospheric wind profile
        du_dtVector = (velocityBallVector(:,iTimeStep)-...
            velocityBallInitial)/dt;%acceleration vector by finite difference
        errorFunction = du_dtVector - ...
            get_acceleration(velocityBallVector(:,iTimeStep),...
            velocityWindAtH,DENSITY_AIR,dragArea,massBall,GRAVITY);
        %difference between LHS and RHS of dynamic equations
        jacobianMatrix = get_jacobian(velocityBallVector(:,iTimeStep),...
            velocityBallInitial,velocityWindAtH,DENSITY_AIR,dragArea,...
            massBall,GRAVITY,dt,errorFunction,jacobianMatrix);
        %numerical Jacobian matrix
        updateVector = solve_equations(jacobianMatrix,errorFunction);
        errorMagnitude = norm(updateVector);
        velocityBallVector(:,iTimeStep) = ...
            velocityBallVector(:,iTimeStep) - updateVector;
    end
    velocityBallInitial = velocityBallVector(:,iTimeStep);
    %updating velocity at time step n to n+1
    positionBallVector(:,iTimeStep) = positionBallVector(:,iTimeStep-1)...
        - velocityBallInitial*dt;%updated position by integrating velocities
    if positionBallVector(2,iTimeStep) < -1.5
        %Simulation stops if the ball reaches 1.5 m below the ground
        positionBallVector(:,iTimeStep+1:end) = [];
        break
    end
end
end

function accelerationVector = get_acceleration(velocityBallVector,...
    velocityWindAtH,DENSITY,dragArea,massBall,GRAVITY)
%This function calculates acceleration vector from horizontal and vertical 
%dynamic equations.
u = velocityBallVector(1);
w = velocityBallVector(2);
velocityRelative = [velocityWindAtH-u;-w];
dragForce = .5*DENSITY*norm(velocityRelative)*dragArea*velocityRelative;
gravityForce = [0;massBall*GRAVITY];
accelerationVector = 1/massBall*(gravityForce+dragForce);
end

function jacobianMatrix = get_jacobian(velocityBallVector,velocityBallInitial...
    ,velocityWindAtH,DENSITY,dragArea,massBall,GRAVITY,dt,errorFunction,jacobianMatrix)
%This function calculates the numerical Jacobian matrix by perturbing each
%of the two velocity parameters.
for i = 1:2
    velocityBallPerturbed = velocityBallVector;
    velocityBallPerturbed(i) = velocityBallVector(i) + 1e-6;
    du_dtVector = (velocityBallPerturbed-velocityBallInitial)/dt;
    errorFunctionPerturbed = du_dtVector - ...
        get_acceleration(velocityBallPerturbed,velocityWindAtH,DENSITY,...
        dragArea,massBall,GRAVITY);
    jacobianMatrix(:,i) = (errorFunctionPerturbed-errorFunction)/1e-6; 
end
end
