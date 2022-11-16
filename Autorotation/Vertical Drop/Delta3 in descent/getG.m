function [loading,G,out] = getG(statesold,states,dt,n,dpsi,psi,model,rotor,physics,airfoil)
M = model.M;
Linv = model.Linv;
Nb = rotor.N_blades;
lambdastates = model.lambda_states;%1 for symmetric Pitt-Peters model

g = physics.g;
rho = physics.rho;
m = rotor.mass;
R = rotor.radius;
dr = 0.1;
r = 0.05:dr:0.95;

omega = states(end-1);
V_desc = states(end);
F = zeros(lambdastates,1);%RHS of dynamic inflow equations
eq34k = zeros(2*Nb,1);%2*N_blades flapping equations
statesstar = (states-statesold)/dpsi;%non-dimensional time derivative
statesdot = (states-statesold)/dt;%dimensional time derivative

lambda0i = states(1); lambdasi = 0*states(2); lambdaci = 0*states(3);
omegadot = statesdot(end-1);
for k = 1:Nb %k: k-th blade; Nb: total # blades
    psik = psi + 2*pi/Nb*(k-1);
    lambdak = lambdasi*r*sin(psik) + lambdaci*r*cos(psik) + lambda0i;% + lambda0;%%
    betaveck = states(lambdastates+(2*k-1:2*k));
    betavecdotk = statesdot(lambdastates+(2*k-1:2*k));
    [eq34k(2*k-1:2*k,1),CTk,CLk,CMk,eq5,Q,alpha,M_flap,Fx1,Fz1] = exactBET(n,...
        lambdak,betaveck,betavecdotk,V_desc,omega,omegadot,r,dr,rotor,physics,airfoil);
    F = F+[CTk;CLk*sin(psik);CMk*cos(psik)];
end
lambda = states(1:lambdastates);
lambdastar = statesstar(1:lambdastates);
lambda_tot = V_desc/omega/R + lambda(1);
l_bar = lambda_tot/sqrt(F(1)/2);
if l_bar > -1 && l_bar < 0.6378
    g_lambda = real((2+l_bar)^-2 - l_bar^2 + (1+l_bar)*(0.109 + 0.217*(l_bar-0.15)^2));
else
    g_lambda = 0;
end
V_T = sqrt(lambda_tot^2 + g_lambda*F(1)/2);
eq123 = M*lambdastar + 1.35*V_T*Linv*lambda - F;
eq6 = -V_desc + statesold(end) - g*dt + F(1)*rho*pi*R^4*(cos(betaveck(2)))^2.*states(end-1)^2/m*dt + 0*0.25*V_desc^2*dt;
out = [alpha,M_flap,Fx1,Fz1];
eq1 = eq123(1);
G = [eq1;eq34k;eq5;eq6];%1,2,3: non-dimensional, 3,4,5: dimensional
loading = [F;Q];
end

function [eq34,CT,CL,CM,eq5,Q,alpha,M_flap,Fy1,Fz1] = exactBET(n,lambda,betavec,betavecdot,...
    V_desc,omega,omegadot,r,dr,rotor,physics,airfoil)
R = rotor.radius;
rho = physics.rho;
mu = physics.mu;
c = rotor.chord;
Ib = rotor.Ib;
delta3 = rotor.delta3;
angle_tw = rotor.twist;

beta = betavec(1);
betadot = betavec(2);
theta = rotor.theta0(n)-beta*tan(delta3)+angle_tw*pi/180*r;
vnormal = (V_desc+lambda*omega*R)*cos(beta)+R*r*betadot*cos(delta3);%blade frame
vtangential = omega*R*r*cos(beta);%blade frame
phi = atan(vnormal./vtangential);
if phi > pi/2; error('error in phi');end
alpha = theta-phi;
V = sqrt(vnormal.^2 + vtangential.^2);
Re = rho*abs(V)*c/mu;
[Cl,Cd] = lookup(alpha,Re,airfoil);
L = .5*rho*V.^2.*c.*Cl;
D = .5*rho*V.^2.*c.*Cd;

% |(z)
% |  /(y)
% | /
% |/_____(x)
Fx1 = -(L.*cos(phi)-D.*sin(phi))*sin(beta);
Fy1 = -(L.*sin(phi)+D.*cos(phi))*cos(beta)+Ib*omega*sin(2*beta)*betadot/R;
Fz1 = (L.*cos(phi)-D.*sin(phi))*cos(beta);

CT = trapz(Fz1)*dr*R/(rho*pi*R^4*omega^2*(cos(beta))^2);%parallel to shaft
CL = trapz(Fz1.*r)*dr/(rho*pi*R^2*omega^2*(cos(beta))^2);%R removed because r non-dimensional
CM = trapz(Fz1.*r)*dr/(rho*pi*R^2*omega^2*(cos(beta))^2);%R removed because r non-dimensional
Q = -trapz(R*r.*L.*sin(phi)*cos(beta))*R*dr-trapz(R*r.*D.*cos(phi)*cos(beta))*R*dr;%driving torque
M_flap = trapz(Fz1.*r*cos(beta) - Fx1.*r*sin(beta))*dr*R^2;%flap up positive

eq3 = Ib*(betavecdot(2) + omega^2*sin(beta)*cos(beta)) - M_flap;% 1st order flap eqn
eq4 = (betavecdot(1) - betavec(2));
eq34 = [eq3;eq4];
eq5 = Ib*omega*sin(2*beta)*betadot + Ib*omegadot*(cos(beta))^2 - Q;% torque eqn
end