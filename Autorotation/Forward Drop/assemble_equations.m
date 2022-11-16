function [loading,G] = assemble_equations(model,rotor,body,physics,airfoil,...
    statesold,states,theta_input,Fx,Fy,U,psi,dpsi,dt)
%This function assembles the system of 9 (or 11) equations
%eq1:       dynamic inflow thrust
%(eq2,3):   dynamic inflow higher harmonics
%eq4 (x2):  blade flap
%eq5 (x2):  betadot
%eq6:       torque
%eq7,8:     horizontal motion
%eq9:       vertical motion and gravity 

% short notations
M = model.M;
Linv = model.Linv;
Nb = rotor.N_blades;
lambdastates = model.lambda_states;%1 for symmetric Pitt-Peters model
betastates = model.beta_states;
R = rotor.radius;
CD_cylinder = physics.CD_cylinder;


% resolving velocities and states
omega = states(lambdastates+betastates+1);
u = states(lambdastates+betastates+2);
v = states(lambdastates+betastates+3);
w = states(lambdastates+betastates+4);
Vh = sqrt((u-U)^2+v^2);
chi = -atan(w/Vh);
eq45k = zeros(2*Nb,1);
statesstar = (states-statesold)/dpsi;
statesdot = (states-statesold)/dt;
lambda0i = states(1); 
lambdasi = states(2)*(lambdastates==3); 
lambdaci = states(3)*(lambdastates==3);
mu = Vh/omega/R;
omegadot = statesdot(lambdastates+betastates+1);

% Initalizations and discretization
eq6 = 0;
F = zeros(3+2,1);
dr = 0.3;
r = 0.05:dr:0.95;%nondimensional radial location

for k = 1:Nb %k: k-th blade; Nb: total # blades
    psik = psi + 2*pi/Nb*(k-1);
    thetak = theta_input(1)+theta_input(2)*cos(psik)+theta_input(3)*sin(psik);
    lambdak = lambdasi*r*sin(psik) + lambdaci*r*cos(psik) + lambda0i;
    betaveck = states(lambdastates+(2*k-1:2*k));
    betavecdotk = statesdot(lambdastates+(2*k-1:2*k));
    [eq45k(2*k-1:2*k,1),CTk,CLk,CMk,eq6k,Q,FX,FY] = exactBET(rotor,...
        physics,airfoil,psik,lambdak,thetak,betaveck,betavecdotk,...
        omega,omegadot,w,Vh,r,dr);%get forces and equations from BET
    % Sum over all blades
    eq6 = eq6+eq6k;%torque equation
    F = F+[CTk;CLk*sin(psik);CMk*cos(psik);FX;FY];
end

% Assembling Dynamic Inflow equations
lambda = states(1:lambdastates);
lambdadot = statesstar(1:lambdastates);
nu = lambda(1);
lambda_h = real(sqrt(F(1)/2));
lambda_tot = w/omega/R + lambda(1);
l_bar = lambda_tot/lambda_h;
if l_bar > -1 && l_bar < 0.6378
    g_lambda = ((2+l_bar)^-2 - l_bar^2 + (1+l_bar)*(0.109 + 0.217*(l_bar-0.15)^2));
    g_prime = real(-2*(2+l_bar)^-3 + 0.049 - 1.696*l_bar + 0.651*l_bar^2);
else
    g_lambda = 0;
    g_prime = 0;
end

mu_bar = mu/lambda_h;
if mu_bar >= 0 && mu_bar < 0.707
    f_mu = 1-2*mu_bar^2;
    f_prime = -4*mu_bar;
else
    f_mu = 0;
    f_prime = 0;
end

V_T = sqrt(lambda_tot^2 + mu^2 + lambda_h^2*g_lambda*f_mu);
V = (lambda_tot*(lambda_tot+nu) + mu^2 + lambda_h^2*(g_lambda+g_prime*nu/lambda_h/2)*f_mu)...
    /(V_T - nu/2*g_lambda*f_mu + nu/4*(lambda_tot/lambda_h*g_prime*f_mu+mu_bar*g_lambda*f_prime));
Vmat = diag([V_T,V,V]);
eq123 = M*lambdadot + Vmat*Linv(pi/2)*lambda - F(1:3);

% Assembling body dynamics equations
eq7 = -u + statesold(lambdastates+betastates+2) + Fx/body.mass*dt - ...
    CD_cylinder*0.5*physics.rho*norm([Vh,w])^2*body.length*body.diameter...
    *cos(chi)^2*(u-U)/Vh/body.mass;% 1e-2*norm(u,v)*u*dt;
eq8 = -v + statesold(lambdastates+betastates+3) + Fy/body.mass*dt - ...
    CD_cylinder*0.5*physics.rho*norm([Vh,w])^2*body.length*body.diameter...
    *cos(chi)^2*v/Vh/body.mass;%1e-2*norm(u,v)*v*dt;
eq9 = -w + statesold(lambdastates+betastates+4) - physics.g*dt + ...
    F(1)*physics.rho*pi*R^4.*omega^2/body.mass*dt + 0*0.25*w^2*dt;

eq1 = eq123(1);
if lambdastates == 3
    G = [eq123;eq45k;eq6;eq7;eq8;eq9];%1,2,3: non-dimensional, others: dimensional
elseif lambdastates == 1
    G = [eq1;eq45k;eq6;eq7;eq8;eq9];%1: non-dimensional, others: dimensional
end
    loading = [F;Q];% CT; CL; CM; FX; FY; Q
end

function [eq45,CT,CL,CM,eq6,Q,FX,FY] = exactBET(rotor,physics,airfoil,...
    psi,lambda,thetak,betavec,betavecdot,omega,omegadot,Vd,Vh,r,dr)
% short notations
R = rotor.radius;
rho = physics.rho;
mu = physics.mu;
c = rotor.chord;
Ib = rotor.Ib;
delta3 = rotor.delta3;

% Calculating blade lift and drag forces
beta = betavec(1);
betadot = betavec(2);
theta = -beta*tand(delta3)+6*pi/180*r*rotor.twist+thetak;%pitch-flap coupling
vnormal = (Vd+lambda*omega*R+Vh*beta)*cos(beta)+R*r*betadot;%blade frame
vtangential = omega*R*r*cos(beta)-Vh*sin(psi);%blade frame
phi = atan(vnormal./vtangential);
if phi > pi/2; error('error in phi');end
alpha = theta-phi;
V = sqrt(vnormal.^2 + vtangential.^2);
Re = rho*V*c/mu;
[Cl,Cd] = lookup(alpha,Re,airfoil);%thin airfoil lookup table
L = .5*rho*V.^2.*c.*Cl;
D = .5*rho*V.^2.*c.*Cd;

% Integrated forces and moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% |(z1) >>>>>> Upwards towards shaft        %%%
% |  /(y) >>>>>> Towards blade leading edge %%%
% | /                                       %%%
% |/_____(x1) >>>>>> Towards blade tip      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fx1 = -(L.*cos(phi)-D.*sin(phi))*sin(beta);
Fy1 = -(L.*sin(phi)+D.*cos(phi));
Fz1 = (L.*cos(phi)-D.*sin(phi))*cos(beta);

FX = trapz(-Fy1*sin(psi)+Fx1*cos(psi))*dr;%ground frame X force
FY = trapz(Fy1*cos(psi)+Fx1*sin(psi))*dr;%ground frame Y force

CT = trapz(Fz1)*dr*R/(rho*pi*R^4*omega^2*(cos(beta))^2);%parallel to shaft
CL = trapz(Fz1.*r)*dr/(rho*pi*R^3*omega^2*(cos(beta))^2);
CM = trapz(Fz1.*r)*dr/(rho*pi*R^3*omega^2*(cos(beta))^2);
Q = -trapz(R*r.*L.*sin(phi)*cos(beta))*R*dr-trapz(R*r.*D.*cos(phi)*cos(beta))*R*dr;%driving torque
M_flap = trapz(Fz1.*r*cos(beta) - Fx1.*r*sin(beta))*dr*R^2;%flap up positive
eq4 = Ib*(betavecdot(2) + omega^2*sin(beta)*cos(beta)) - M_flap;%flap equation 1
eq5 = betavecdot(1) - betavec(2);%flap equation 2
eq45 = [eq4;eq5];
eq6 = Ib*omega*sin(2*beta)*betadot + Ib*omegadot*cos(beta)^2 -  Q;%torque equation
end