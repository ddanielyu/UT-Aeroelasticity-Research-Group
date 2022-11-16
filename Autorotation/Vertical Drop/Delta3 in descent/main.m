clearvars;
% close all;
clc;
saving = false;
plotting = true;
tic

%% Read Airfoil tables
airfoil.type = 'arc';%'sym' or 'arc'
load('aeroTable');%circular arc
airfoil.arc.alpha = aeroTable(:,1);%circular arc
airfoil.arc.Cl = aeroTable(:,2);%circular arc
airfoil.arc.Cd = aeroTable(:,3);%circular arc

load('symmetric_airfoil.mat')%symmetric
airfoil.sym.al = struct('al1',al1,  'al2',al2, 'al4',al4, 'al8',al8,... 
    'al16',al16, 'al36',al36, 'al70',al70, 'al100',al100, 'al200',al200,...
    'al500',al500, 'al1000',al1000);
airfoil.sym.cl = struct('cl1',cl1,  'cl2',cl2, 'cl4',cl4, 'cl8',cl8,... 
    'cl16',cl16, 'cl36',cl36, 'cl70',cl70, 'cl100',cl100, 'cl200',cl200,...
    'cl500',cl500, 'cl1000',cl1000);
airfoil.sym.cd = struct('cd1',cd1,  'cd2',cd2, 'cd4',cd4, 'cd8',cd8,... 
    'cd16',cd16, 'cd36',cd36, 'cd70',cd70, 'cd100',cd100, 'cd200',cl200,...
    'cd500',cd500, 'cd1000',cd1000);


%% Simulation parameters
tsteps = 201;
logsteps = 70;%log increments up to 50 steps, then constant
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
rotor.radius = .22;
rotor.chord = 0.025;
rotor.mass = 0.15;
rotor.Ib = 0.000151;
rotor.delta3 = 30*pi/180;
rotor.twist = 6;
rotor.theta0 = 0*ones(tsteps,1);

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
    n
    psi = psi+dpsi;
    dpsi = dt(n)*statesold(end-1);
    iter = 0;
    err = 1;
    while err > hyperparameter.tolerance && iter < hyperparameter.max_iterations
        iter = iter+1;
        [loading,G,out] = getG(statesold,states,dt(n),n,dpsi,psi,model,rotor,physics,airfoil);%%%
        Jac = getJac(statesold,states,dt(n),n,dpsi,psi,G,model,rotor,physics,airfoil,hyperparameter);%%%
        
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
    err;
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
end
toc

%%
function Jac = getJac(statesold,states,dt,n,dpsi,psi,G,model,rotor,physics,airfoil,hyperparameter)
epsilon = hyperparameter.epsilon;
S = model.total_states;
Jac = zeros(S,S);
for i = 1:S
    statespert = states;
    statespert(i) = states(i)+epsilon;
    [~,G_pert,~] = getG(statesold,statespert,dt,n,dpsi,psi,model,rotor,physics,airfoil);
    Jac(:,i) = (G_pert-G)/epsilon;
    
end
end