function [Cl,Cd] = lookup(alpha_rad,Re,airfoil)
alpha_fluent = airfoil.alpha;
cl_fluent = airfoil.Cl;
cd_fluent = airfoil.Cd;

alpha_deg = alpha_rad*180/pi;
Cl=1.1*interp1(alpha_fluent,cl_fluent,alpha_deg,'linear');
Cd=interp1(alpha_fluent,cd_fluent,alpha_deg,'linear');
end