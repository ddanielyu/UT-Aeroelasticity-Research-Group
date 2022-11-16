function eq789 = global_equations()
%---------------------------%
% North:x | East:y | Down:z %
%---------------------------%

v_vec = [u;v;w];
v_old_vec = [];
v_dot_vec = (v_vec-v_old_vec)/dt;
v_wind_vec = [U;V;0];
v_rel_vec = v_wind_vec-v_vec;
v_rel = norm(v_rel_vec);
area_ref = body.length*body.diameter;

alpha = atan(w/norm(v_rel_vec(1:2)));
CD = 1.1*sin(alpha).^3+0.2*cos(alpha).^3;
CL = 1.1*sin(alpha).^2.*cos(alpha)+0.2*cos(alpha).^2.*sin(alpha);

force_gravity_vec = [0;0;rotor.mass*physics.g];
force_rotor_vec = [Fx;Fy;Fz];
force_body_vec = [];
eq9 = -w + statesold(lambdastates+betastates+4) - physics.g*dt + ...
    F(1)*physics.rho*pi*R^4*(cos(betaveck(2)))^2.*omega^2/rotor.mass*dt + 0*0.25*w^2*dt;


eq789 = v_dot_vec - force_gravity_vec - force_rotor_vec - force_body_vec;
% [eq7;eq8;eq9];

end
% 
% function [Cd,Cl] = cylinder_coefficients(alpha, Re, L, D)
% %% Lift coefficient curve fit
% l1 = 1.9276;
% l2 = -1.9124;
% l3 = 0.4270;
% l4 = 0.3870;
% l5 = -0.4420;
% l6 = 0.5050;
% l7 = 0.0330;
% l8 = -0.1366;
% l9 = 0.3912;
% l10 = 0.0515;
% l11 = -0.0588;
% l12 = 0.0382;
% 
% d1 = l1 + l2*exp(-l3*L/D);
% d2 = l4 + l5*exp(-l6*L/D);
% d3 = l7 + l8*exp(-l9*L/D);
% d4 = l10 + l11*exp(-l12*L/D);
% 
% A2 = d1/Re + d2;
% A4 = A2*(d3*log(Re) + d4);
% 
% Cl = A2*sin(2*alpha) + A4*sin(4*alpha);
% 
% %% Drag coefficient curve fit
% % clearvars -except alpha Re L D
% g1 = -0.0414;
% g2 = 0.07010;
% g3 = 0.9342;
% g4 = -0.2377;
% g5 = 0.2685;
% g6 = 0.1697;
% g7 = -0.04070;
% g8 = 0.0790;
% g9 = 0.8194;
% g10 = 0.7609;
% g11 = 0.2577;
% g12 = 0.1583;
% 
% e1 = 18.7565;
% e2 = 10.8586;
% e3 = -0.3233;
% e4 = -0.6024;
% 
% b1 = g1 + g2*exp(-g3*L/D);
% b2 = g4 + g5*exp(-g6*L/D);
% b3 = g7 + g8*exp(-g9*L/D);
% b4 = g10 + g11*exp(-g12*L/D);
% 
% k1 = e1*D/L + e2;
% k2 = e3*D/L + e4;
% 
% A2 = b1*log(Re) + b2;
% A0 = b3*log(Re) + b4;
% 
% Cd_perp = k1*Re^k2;
% Cd = Cd_perp*(A2*cos(2*alpha) + A0);
% 
% end