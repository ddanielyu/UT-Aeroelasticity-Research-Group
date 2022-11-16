function errorUpdateTruncated = solve_equations(jacobianMatrix,rhsVector,...
    varargin)
%solves system of linearized equations
%if indices specified: drop them in jacobianMatrix, rhsVector and errorUpdate
%if rhsVector below tolerance: drop corresponding indices
tolerance = 1e-12;
lengthG = length(rhsVector);
errorUpdateTruncated = zeros(lengthG,1);
iSmall = find(abs(rhsVector) < tolerance);
iDropped = union(iSmall,cell2mat(varargin));

rhsVector(iDropped) = [];
jacobianMatrix(iDropped,:) = [];
jacobianMatrix(:,iDropped) = [];
errorUpdate = jacobianMatrix\rhsVector;
errorUpdateTruncated(setdiff(1:lengthG,iDropped)) = errorUpdate;

end