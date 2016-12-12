function [r2,Niter] = minimize_r(atoms,bonds,r1,DBCIndices,tDBCValues,...
    FreeIndices,R0,z1,C,TOL_r)

% Solve for x, a vector of free degrees of freedom, using standard Newton algorithm
x = r1(FreeIndices);
eps_r = 1+TOL_r;
Niter = 0;
while eps_r > TOL_r
    Niter = Niter+1;
    [r2,f_r,K_r] = grad_hess(atoms,bonds,x,DBCIndices,tDBCValues,...
        FreeIndices,R0,z1);
    
    % Introduce constraints, use primal-dual formulation, assembly extended quantities
    EK_r = [K_r,C(:,FreeIndices)'
        C(:,FreeIndices),sparse(size(C,1),size(C,1))];
    Ef_r = [f_r;C*r2];
    
    % Solve the system
    du = -EK_r\Ef_r;
    dx = du(1:end-size(C,1));
    lambda = du(end-size(C,1)+1:end);
    x = x+dx;
    
    % Update the error
    eps_r = norm(dx)/norm(x-R0(FreeIndices))+...
        norm(f_r+C(:,FreeIndices)'*lambda);
end

% Reconstruct the converged r-vector from the iterative variable x
r2(DBCIndices) = R0(DBCIndices)+tDBCValues;
r2(FreeIndices) = x;

end
