function er = eq_solve(Jac,G,varargin)
%if indices specified, drop G, Jac and er
%if G below tolerance, drop G, Jac and er
tolerance = 1e-12;
length_G = length(G);
er = zeros(length_G,1);
ind_small = find((abs(G) < tolerance)+ isnan(G));
ind_dropped = union(ind_small,cell2mat(varargin));
% G(ind_dropped)
G(ind_dropped) = [];
Jac(ind_dropped,:) = [];
Jac(:,ind_dropped) = [];
ermod = Jac\G;
er(setdiff(1:length_G,ind_dropped)) = ermod;

end