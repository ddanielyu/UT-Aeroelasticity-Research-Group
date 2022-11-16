


if ~exist('main_flag','var')
    clc;
    clearvars;
    plotting = true;
    tic
end
saving = false;
%% Read Airfoil tables
load('aeroTable')
airfoil.alpha = aeroTable(:,1);
airfoil.Cl = aeroTable(:,2);
airfoil.Cd = aeroTable(:,3);

%% Simulation parameters
tsteps = 201;
logsteps = 50;%log increments up to 50 steps, then constant
dt = 0.01*pi*[(logspace(-2,0,logsteps)),ones(1,tsteps-logsteps)];
time_vec = cumtrapz(dt);

hyperparameter.epsilon = 1e-6;
hyperparameter.tolerance = 1e-12;
hyperparameter.max_iterations = 50;

%% Physical constants
physics.rho = 1.225;
physics.g = 9.8;
physics.mu = 1.73e-5;

%% Rotor properties
rotor.N_blades = 2;
rotor.radius =  convlength(2.25,'ft','m');%convlength(23.6/2,'in','m');%
if ~exist('main_flag','var')
    rotor.chord = convlength(1.126,'in','m');%0.025/.22*rotor.radius;%
    rotor.delta3 = 30*pi/180;
    rotor.twist = 0;
    rotor.theta0 = 0;
else
    rotor.chord = rotor.radius/AR;
    rotor.delta3 = delta3;
    rotor.twist = theta_tw;
    rotor.theta0 = theta0;
    plotting = false;
end
rotor.mass = convmass(5.5,'lbm','kg');%.675;%
rotor.Ib = 0.000151/.22^5*rotor.radius^5;%5e-3*rotor.radius^2/3;%0.000151/.22^5*rotor.radius^5;%

%% Model properties
model.M = [128/75/pi,0,0;0,-256/945/pi,0;0,0,-256/945/pi];%
L = [1/2,0,0;0,-2,0;0,0,-2];
model.Linv = inv(L);
model.lambda_states = 1;
model.total_states = model.lambda_states+2*rotor.N_blades+2;
%lambda0,beta1,betadot1,beta2,betadot2,omega,descent_speed

%% Initialize
beta0 = 0;%initial blade flapping angle
statesold = [zeros(model.total_states-2*rotor.N_blades-2,1);repmat([beta0;0],rotor.N_blades,1);2*pi;-.1];%1 Hz starting RPM
F = zeros(model.lambda_states+1,tsteps);%Dimensional thrust and torque
psi = 0;
dpsi = 0;

%% Assume
states = statesold;%initial assumption for the next time step
states_theoretical = statesold*ones(1,tsteps);
flag = 0;%0: unlocked; 1: locked
output = zeros(31,length(states_theoretical));
%% Update/March and Iterate
for n = 2:tsteps%t=0 is initial known steady state%%%
    n;
    psi = psi+dpsi;
    dpsi = dt(n)*statesold(end-1);
    iter = 0;
    err = 1;
    while err > hyperparameter.tolerance && iter < hyperparameter.max_iterations
        iter = iter+1;
        [loading,G,out] = getG(statesold,states,dt(n),dpsi,psi,model,rotor,physics,airfoil);%%%
        Jac = getJac(statesold,states,dt(n),dpsi,psi,G,model,rotor,physics,airfoil,hyperparameter);%%%
        
        if flag == 0 && statesold(2) >= pi/4 % unlock -> lock
            %             flag
            Jacmod = Jac([1,6:end],[1,6:end]);
            Gmod = G([1,6:end]);
            er = Jacmod\Gmod;
            er = [er(1:model.lambda_states);0;0;0;0;er(end-1:end)];
            states([3,5]) = [0,0];
            if G(2) > hyperparameter.tolerance
                flag = 1; %lock -> unlock
            end
        else
            flag = 0; %unlock
            %             flag = 'inflow';
            er = Jac\G;
        end
        err = norm(er);
        states = states - er;
        
    end
    %     err
    output(:,n) = out;%All blade distributed parameters from BET go here
    states_theoretical(:,n) = states;%States output
    statesold = states; %updating time step n to n+1
    F(:,n) = loading([1;end]);%Loads output
end
if saving
    save('drop_simulation.mat') %#ok<*UNRCH>
end
if plotting
    plotfun;
    toc
end

%%
function Jac = getJac(statesold,states,dt,dpsi,psi,G,model,rotor,physics,airfoil,hyperparameter)
epsilon = hyperparameter.epsilon;
S = model.total_states;
Jac = zeros(S,S);
for i = 1:S
    statespert = states;
    statespert(i) = states(i)+epsilon;
    [~,G_pert,~] = getG(statesold,statespert,dt,dpsi,psi,model,rotor,physics,airfoil);
    Jac(:,i) = (G_pert-G)/epsilon;
    
end
end