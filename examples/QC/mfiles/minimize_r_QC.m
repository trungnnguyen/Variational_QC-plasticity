function [r2,r,Niter] = minimize_r_QC(atoms,samplingatoms,bonds,r1,...
    DBCIndicesQC,tDBCValuesQC,FreeIndicesQC,R0QC,z1,C,Phi,TOL_r)

% Solve for x - a vector of free degrees of freedom, using standard Newton algorithm
x = r1(FreeIndicesQC);
eps_r = 1+TOL_r;
Niter = 0;
while eps_r > TOL_r
    Niter = Niter+1;
    [r2,G,H] = grad_hess_QC(atoms,samplingatoms,bonds,x,DBCIndicesQC,...
        tDBCValuesQC,FreeIndicesQC,R0QC,z1,Phi);
    
    % Introduce constraints, use primal-dual formulation, assembly extended quatities
    EH = [H,C(:,FreeIndicesQC)'
        C(:,FreeIndicesQC),sparse(size(C,1),size(C,1))];
    EG = [G;C*r2];
    
    % Solve the system
    du = -EH\EG;
    dx = du(1:end-size(C,1));
    lambda = du(end-size(C,1)+1:end);
    x = x+dx;
    
    % Update the error
    eps_r = norm(dx)/norm(x-R0QC(FreeIndicesQC))+...
        norm(G+C(:,FreeIndicesQC)'*lambda);
end

% Reconstruct converged r2- and r-vectors from x
r2 = zeros(size(R0QC,1),1);
r2(DBCIndicesQC) = R0QC(DBCIndicesQC)+tDBCValuesQC;
r2(FreeIndicesQC) = x; % QC displacement vector reconstruction
r = Phi*r2; % full r reconstruction

end
