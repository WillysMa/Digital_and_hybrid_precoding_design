function [y, cost] = MO_Frf_LS(FRF, Fbb,Fopt)
[Nt, NRF] = size(FRF);

% importmanopt()
manifold = complexcirclefactory(Nt, NRF);
problem.M = manifold;

problem.cost = @(X)  norm(Fopt-X*Fbb,'fro')^2;
problem.egrad = @(X) 2*X*Fbb*Fbb'-2*Fopt*Fbb';

% checkgradient(problem);
warning('off', 'manopt:getHessian:approx');

[x,cost,info,options] = conjugategradient(problem,FRF);
% [x,cost,info,options] = trustregions(problem, FRF(:));
y = reshape(x,Nt,NRF);

ccc=1;
end
