clearvars;
% close all;
clc;
saving = false;
plotting = 'none';%'all', 'none' or 'OutputsOnly'
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
psisteps = 4001;
logsteps = 50;%log increments up to 50 steps, then constant
dpsi = 0.1*pi*[(logspace(-2,0,logsteps)),ones(1,psisteps-logsteps)];
psi_vec = cumsum(dpsi);
time_vec = zeros(size(psi_vec));
[~,idx] = findpeaks(-diff(mod(psi_vec,2*pi)),'MinPeakHeight',3);%finds instances of complete rotations
idx = [0,idx,psisteps];
n_rot = 2;

hyperparameter.epsilon = 1e-6;
hyperparameter.tolerance = 1e-12;
hyperparameter.max_iterations = 50;

%% Physical constants
physics.rho = 1.225;
physics.mu = 1.73e-5;
physics.g = 9.8;
physics.CD_cylinder = 0.4;%Hoerner page 3-11

%% Rotor properties
rotor.N_blades = 2;
rotor.radius = .321;%convlength(23.6/2,'in','m');%.22;%24.5
rotor.chord = .030;%convlength(1.126,'in','m');%0.0286;%1.1875
rotor.Ib = 6.5*1e-3*rotor.radius^2/3;%0.000151;%4.5*%16.5
rotor.delta3 = 30*pi/180;
rotor.twist = 0;%0=> untwisted blades, 1=> twisted blades 6 deg root to tip

%% Body properties
body.mass = convmass(1.5,'lbm','kg');
body.length = 0.5;%m
body.diameter = 0.1;%m


%% Control
%steady
control.theta0bar = -0*pi/180;
control.theta1cbar = 0;
control.theta1sbar = 3*pi/180*ones(1,length(psi_vec));
control.beta0 = pi/4;%initial flap angle or precone

%ramp
logsteps = 120;
ramplength = 40;
control.theta0 = 0*[zeros(1,logsteps),linspace(0,15*pi/180,ramplength),...
    15*pi/180*ones(1,psisteps-logsteps-ramplength)];
control.theta1c = [zeros(1,logsteps),linspace(0,1*pi/180,ramplength),...
    1*pi/180*ones(1,psisteps-logsteps-ramplength)];
control.theta1s = 0*[zeros(1,logsteps),linspace(0,1*pi/180,ramplength),...
    1*pi/180*ones(1,psisteps-logsteps-ramplength)];

%Wind
U_init = 5;%mph
buffer = 200;%No wind for buffer time steps
U_10 = [zeros(buffer,1);markov_chain(psisteps-buffer,U_init+1)]*0.44704;%m/s @10m
U_vec = U_10;%initialize
h1 = 10;%reference data available at 10 m
h0 = 0.4;%https://websites.pmc.ucsc.edu/~jnoble/wind/extrap/
altitude = 50;%drop altitude

%% Model properties
model.M = [128/75/pi,0,0;0,-256/945/pi,0;0,0,-256/945/pi];
model.Linv = @(chi) inv([1/2,0,15*pi/64*sqrt((1-sin(chi))/(1+sin(chi)));...
    0,4/(1+sin(chi)),0;...
    15*pi/64*sqrt((1-sin(chi))/(1+sin(chi))),0,4*sin(chi)/(1+sin(chi))]);%inv(L);
model.lambda_states = 1;
model.beta_states = 2*rotor.N_blades;
model.velocity_states = 3;
model.total_states = model.lambda_states+model.beta_states+1+model.velocity_states;
%lambda0, (lambda1s,lambda1c,) beta1,betadot1,beta2,betadot2,omega,u,v,w

%% Initialize
statesold = [zeros(model.lambda_states,1);...
    repmat([control.beta0;0],rotor.N_blades,1);2*pi;0;0;-0.1];%1 Hz starting RPM
F = zeros(6,psisteps);%Dimensional thrust and torques
Jac = zeros(model.total_states);
Jac_save = zeros(model.total_states,model.total_states,psisteps);
Fx = zeros(psisteps,1);
Fy = Fx;
beta0 = Fx;
beta1s = Fx;
beta1c = Fx;

%% Assume
states = statesold;%initial assumption for the next time step
states_theoretical = statesold*ones(1,psisteps);
flag = 0;%0: unlocked; 1: locked
output = zeros(31,length(states_theoretical));

%% Update/March and Iterate
for n = 2:psisteps%t=0 is initial known steady state
    n
    %     control.Vd = Vd(n);
    psin = psi_vec(n);
    dpsin = dpsi(n);
    dt = dpsin/statesold(model.lambda_states+model.beta_states+1);
    time_vec(n) = time_vec(n-1) + dt;
    theta_input = [control.theta0bar+control.theta0(n),...
        control.theta1cbar+control.theta1c(n),...
        control.theta1sbar+control.theta1s(n)];
    U_vec(n) = U_10(n)*log(altitude/h0)/log(h1/h0);
    
    iter = 0;
    err = 1;
    while err > hyperparameter.tolerance && iter < 30
        iter = iter+1;
        
        [loading,G] = assemble_equations(model,rotor,body,physics,airfoil,statesold,...
            states,theta_input,Fx(n-1),Fy(n-1),U_vec(n),psin,dpsin,dt);
        
        Jac = getJac(model,rotor,body,physics,airfoil,hyperparameter,...
            statesold,states,theta_input,Fx(n-1),Fy(n-1),U_vec(n),psin,dpsin,dt,G,Jac);
        tempG = G(model.lambda_states+1);
        tempbeta = statesold(model.lambda_states+1)*180/pi;
        flag;
        if flag == 0 && (statesold(model.lambda_states+1) >= 0.95*pi/4 ||...
                statesold(model.lambda_states+3) >= 0.95*pi/4)
            states(model.lambda_states+[1,3]) = [pi/4,pi/4];
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
    %     stak(:,n) = [out,G(4)];
    states([model.lambda_states+1,model.lambda_states+3]) = min(states([model.lambda_states+1,model.lambda_states+3]),pi/4);
    states_theoretical(:,n) = states;%output
    statesold = states; %updating time step n to n+1
    F(:,n) = loading;%output
    
    %Integrate over rotations
    if n == idx(n_rot)
        Fx((idx(n_rot)+1):end) = trapz(F(4,(idx(n_rot-1)+1):idx(n_rot))...
            .*dpsi((idx(n_rot-1)+1):idx(n_rot)));
        Fy((idx(n_rot)+1):end) = trapz(F(5,(idx(n_rot-1)+1):idx(n_rot))...
            .*dpsi((idx(n_rot-1)+1):idx(n_rot)));
        beta0((idx(n_rot)+1):end) = 1/2/pi*trapz(states_theoretical...
            (model.lambda_states+1,(idx(n_rot-1)+1):idx(n_rot)).*dpsi((idx(n_rot-1)+1):idx(n_rot)));
        beta1s((idx(n_rot)+1):end) = 1/pi*trapz(states_theoretical(model.lambda_states+1,...
            (idx(n_rot-1)+1):idx(n_rot)).*sin(psi_vec((idx(n_rot-1)+1):idx(n_rot))).*dpsi((idx(n_rot-1)+1):idx(n_rot)));
        beta1c((idx(n_rot)+1):end) = 1/pi*trapz(states_theoretical(model.lambda_states+1,...
            (idx(n_rot-1)+1):idx(n_rot)).*cos(psi_vec((idx(n_rot-1)+1):idx(n_rot))).*dpsi((idx(n_rot-1)+1):idx(n_rot)));
        
        n_rot = n_rot+1;
    end
    ht = 100+cumtrapz(states_theoretical(end,:).*dt);
    altitude = ht(end);
    
    
    
end
figure()
plot3(cumtrapz(time_vec,states_theoretical(end-2,:)),cumtrapz(time_vec,states_theoretical(end-1,:)),ht-ht(end));axis equal
title('Trajectory of Descent')
box on

%save('1600tsteps.mat')
% plotfun;
%%
toc

function Jac = getJac(model,rotor,body,physics,airfoil,hyperparameter,...
    statesold,states,theta_input,Fx,Fy,U,psi,dpsi,dt,G,Jac)
for i = 1:model.total_states
    statespert = states;
    statespert(i) = states(i)+hyperparameter.epsilon;
    [~,G_pert] = assemble_equations(model,rotor,body,physics,airfoil,statesold,...
        statespert,theta_input,Fx,Fy,U,psi,dpsi,dt);
    Jac(:,i) = (G_pert-G)/hyperparameter.epsilon;
    
end
end
%}