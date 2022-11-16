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
psisteps = 201;
logsteps = 50;%log increments up to 100 steps, then constant
xvec = 1*[(logspace(-2,0,logsteps)),ones(1,psisteps-logsteps)];
xvec = [xvec(1),xvec(2:end)];
dpsi = xvec*pi;
psi_vec = cumsum(dpsi);
time_vec = zeros(size(psi_vec));

hyperparameter.epsilon = 1e-6;
hyperparameter.tolerance = 1e-12;
hyperparameter.max_iterations = 50;

%% Physical constants
physics.rho = 1.225;
physics.mu = 1.73e-5;

%% Rotor properties
rotor.N_blades = 2;
rotor.radius = convlength(23.6/2,'in','m');%.22;%24.5
rotor.chord = convlength(1.126,'in','m');%0.0286;%1.1875
rotor.Ib = 6.5*1e-3*rotor.radius^2/3;%0.000151;%4.5*%16.5
rotor.delta3 = 30*pi/180;
rotor.twist = 0;%0=> untwisted blades, 1=> twisted blades 6 deg root to tip

%% Control
%steady
control.chi = pi/2;%degrees rotor axis orientation wrt Vh
Vtunnel = 4;%m/s
control.Vd = -Vtunnel*sin(control.chi);%known dimensional descent velocity m/s
control.Vh = Vtunnel*cos(control.chi);%known dimensional horizontal velocity m/s
control.theta0bar = 0*pi/180;
control.theta1cbar = 0;
control.theta1sbar = 0;
control.beta0 = pi/4;%initial flap angle or precone

%ramp
startstep = 0;
endstep = 60;
stepsize = 6;
control.theta0 = 0*[linspace(startstep,endstep,stepsize),endstep*ones(1,psisteps-stepsize)];
control.theta1c = 0*3*pi/180*ones(1,length(psi_vec));%zeros(1,length(psi));
control.theta1s = 0*3*pi/180*ones(1,length(psi_vec));

%% Model properties
model.M = [128/75/pi,0,0;0,-256/945/pi,0;0,0,-256/945/pi];%
L = [1/2,0,15*pi/64*sqrt((1-sin(control.chi))/(1+sin(control.chi)));...
    0,4/(1+sin(control.chi)),0;...
    15*pi/64*sqrt((1-sin(control.chi))/(1+sin(control.chi))),0,4*sin(control.chi)/(1+sin(control.chi))];
model.Linv = inv(L);
model.lambda_states = 1;
model.beta_states = 2*rotor.N_blades;
model.total_states = model.lambda_states+model.beta_states+1;
%lambda0,beta1,betadot1,beta2,betadot2,omega

%% Initialize
statesold = [zeros(model.lambda_states,1);repmat([control.beta0;0],rotor.N_blades,1);2*pi];%1 Hz starting RPM
F = zeros(6,psisteps);%Dimensional thrust and torques
Jac = zeros(model.total_states,model.total_states);
Jac_save = zeros(model.total_states,model.total_states,psisteps);

%% Assume
states = statesold;%initial assumption for the next time step
states_theoretical = statesold*ones(1,psisteps);
flag = 0;%0: unlocked; 1: locked
output = zeros(31,length(states_theoretical));

%% Update/March and Iterate
for n = 2:psisteps%t=0 is initial known steady state
    n
    psin = psi_vec(n);
    dpsin = dpsi(n);
    dt = dpsin/statesold(model.lambda_states+model.beta_states+1);
    time_vec(n) = time_vec(n-1) + dt;
    theta_input = [control.theta0bar+control.theta0(n),...
        control.theta1cbar+control.theta1c(n),...
        control.theta1sbar+control.theta1s(n)];
    iter = 0;
    err = 1;
    while err > hyperparameter.tolerance && iter < 30
        iter = iter+1;
        [loading,G,out] = getG(model,rotor,control,physics,airfoil,statesold,...
            states,theta_input,psin,dpsin,dt);

        Jac = getJac(model,rotor,control,physics,airfoil,hyperparameter,statesold,states,theta_input,psin,dpsin,dt,G,Jac);
        tempG = G(model.lambda_states+1);
        tempbeta = statesold(model.lambda_states+1)*180/pi;
        flag;
        if flag == 0 && (statesold(model.lambda_states+1) >= 0.95*pi/4 ||...
                statesold(model.lambda_states+3) >= 0.95*pi/4)
            states([2,4]) = [pi/4,pi/4];
            er = eq_solve(Jac,G,model.lambda_states+[1,2,3,4]);
            states([model.lambda_states+2,model.lambda_states+4]) = [0,0];
            if G(model.lambda_states+1) > hyperparameter.tolerance
                flag = 1 %blade flaps back to unlock position
            end
        else
            er = eq_solve(Jac,G);
        end
        Jac_save(:,:,n) = Jac;
        err = norm(er);
        states = states - er;
    end
    iter;
    if isnan(err)
        return
    end
    stak(n) = out;
    states([model.lambda_states+1,model.lambda_states+3]) = min(states([model.lambda_states+1,model.lambda_states+3]),pi/4);
    states_theoretical(:,n) = states;%output
    statesold = states; %updating time step n to n+1
    F(:,n) = loading;%output
end
plotfun;
%%
toc

function Jac = getJac(model,rotor,control,physics,airfoil,hyperparameter,statesold,states,theta_input,psi,dpsi,dt,G,Jac)
for i = 1:model.total_states
    statespert = states;
    statespert(i) = states(i)+hyperparameter.epsilon;
    [~,G_pert,~] = getG(model,rotor,control,physics,airfoil,statesold,statespert,theta_input,psi,dpsi,dt);
    Jac(:,i) = (G_pert-G)/hyperparameter.epsilon;
    
end
end
%}