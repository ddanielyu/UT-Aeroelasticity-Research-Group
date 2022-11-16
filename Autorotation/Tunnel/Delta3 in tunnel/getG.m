function [loading,G,out] = getG(model,rotor,control,physics,airfoil,statesold...
    ,states,theta_input,psi,dpsi,dt)%,thetadot,omega)
M = model.M;
Linv = model.Linv;
Nb = rotor.N_blades;
lambdastates = model.lambda_states;%1 for symmetric Pitt-Peters model
betastates = model.beta_states;
R = rotor.radius;
dr = 0.1;
r = 0.05:dr:0.95;

omega = states(lambdastates+betastates+1);
F = zeros(3+2,1);
eq45k = zeros(2*Nb,1);
statesstar = (states-statesold)/dpsi;
statesdot = (states-statesold)/dt;
lambda0i = states(1); 
lambdasi = states(2)*(lambdastates==3); 
lambdaci = states(3)*(lambdastates==3);
mu = control.Vh/omega/R;
omegadot = statesdot(lambdastates+betastates+1);
eq6 = 0;
for k = 1:Nb %k: k-th blade; Nb: total # blades
    psik = psi + 2*pi/Nb*(k-1);
    thetak = theta_input(1)+theta_input(2)*cos(psik)+theta_input(3)*cos(psik);
    lambdak = lambdasi*r*sin(psik) + lambdaci*r*cos(psik) + lambda0i;
    betaveck = states(lambdastates+(2*k-1:2*k));
    betavecdotk = statesdot(lambdastates+(2*k-1:2*k));
    [eq45k(2*k-1:2*k,1),CTk,CLk,CMk,eq6k,Q,out] = exactBET(rotor,control,physics,airfoil,psik,lambdak,...
        thetak,betaveck,betavecdotk,omega,omegadot,r,dr);
    eq6 = eq6+eq6k;
    F = F+[CTk;CLk*sin(psik);CMk*cos(psik);out.FX;out.FY];%*cos(psik)];
end
lambda = states(1:lambdastates);
lambdadot = statesstar(1:lambdastates);
nu = lambda(1);
lambda_h = real(sqrt(F(1)/2));
lambda_tot = control.Vd/omega/R + lambda(1);
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
k = 2;%momentum correction factor
eq123 = M*lambdadot + k*Vmat*Linv*lambda - F(1:3);

eq1 = eq123(1);
if lambdastates == 3
    G = [eq123;eq45k;eq6];%1,2,3: non-dimensional, 4,5,6: dimensional
elseif lambdastates == 1
    G = [eq1;eq45k;eq6];%1: non-dimensional, 4,5,6: dimensional
end
    loading = [F;Q];
end

function [eq45,CT,CL,CM,eq6,Q,out] = exactBET(rotor,control,physics,airfoil,psi,lambda,thetak,betavec,betavecdot,...
    omega,omegadot,r,dr)
R = rotor.radius;
rho = physics.rho;
mu = physics.mu;
c = rotor.chord;
Ib = rotor.Ib;
delta3 = rotor.delta3;
Vd = control.Vd;
Vh = control.Vh;
%returns non-dimensional values
% mu = Vh;%/omega/R;
beta = betavec(1);
betadot = betavec(2);
theta = -beta*tand(delta3)+thetak;%theta1s*sin(psi)+theta1c*cos(psi);
vnormal = (Vd+lambda*omega*R+0*Vh*beta)*cos(beta)+R*r*betadot;%blade frame
vtangential = omega*R*r*cos(beta)-0*Vh*sin(psi);%blade frame
phi = atan(vnormal./vtangential);
if phi > pi/2; error('error in phi');end;
alpha = theta+6*pi/180*r*rotor.twist-phi;
V = sqrt(vnormal.^2 + vtangential.^2);
Re = rho*V*c/mu;
[Cl,Cd] = lookup(alpha,Re,airfoil);
L = .5*rho*V.^2.*c.*Cl;
D = .5*rho*V.^2.*c.*Cd;

% |(z)
% |  /(y)
% | /
% |/_____(x)
Fx1 = -(L.*cos(phi)-D.*sin(phi))*sin(beta);
Fy1 = -(L.*sin(phi)+D.*cos(phi));
Fz1 = (L.*cos(phi)-D.*sin(phi))*cos(beta);

FX = trapz(Fy1*cos(psi)+Fx1*sin(psi))*dr;
FY = trapz(Fy1*sin(psi)-Fx1*cos(psi))*dr;

CT = trapz(Fz1)*dr*R/(rho*pi*R^4*omega^2*(cos(beta))^2);%parallel to shaft
CL = trapz(Fz1.*r)*dr/(rho*pi*R^3*omega^2*(cos(beta))^2);%R removed because r non-dimensional
CM = trapz(Fz1.*r)*dr/(rho*pi*R^3*omega^2*(cos(beta))^2);%R removed because r non-dimensional
Q = -trapz(R*r.*L.*sin(phi)*cos(beta))*R*dr-trapz(R*r.*D.*cos(phi)*cos(beta))*R*dr;%driving torque
M_flap = trapz(Fz1.*r*cos(beta) - Fx1.*r*sin(beta))*dr*R^2;%flap up positive
eq4 = Ib*(betavecdot(2) + omega^2*sin(beta)*cos(beta)) - M_flap;
eq5 = betavecdot(1) - betavec(2);
eq45 = [eq4;eq5];
eq6 = Ib*omega*sin(2*beta)*betadot + Ib*omegadot*cos(beta)^2 -  Q;
out.FX = FX;
out.FY = FY;
out.Fy1 = Fy1;
out.Fz1 = Fz1;
out.alpha = alpha;
out.Re = Re;
end