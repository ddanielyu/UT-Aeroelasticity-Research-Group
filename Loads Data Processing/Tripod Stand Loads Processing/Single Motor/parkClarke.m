function Q = parkClarke(A,B,C)
%This function calculates the quadrature (Q) current from A,
%B, and C phase currents by assuming all current produced from controller
%is IQ
%{
Inputs: 
A,B,C -- phase currents (column vectors)

Outputs:
Q -- Quadrature current (column vector)
%}

for i = 1:length(A)
    Q(i) = sqrt(2/3)*norm([A(i) B(i) C(i)]);
end
Q = Q';

end

