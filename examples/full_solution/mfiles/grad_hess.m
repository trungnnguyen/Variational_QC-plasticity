function [r2,f_r,K_r] = grad_hess(atoms,bonds,x,DBCIndices,tDBCValues,...
    FreeIndices,R0,z)

% Reconstruct r2
r2 = zeros(2*length(atoms),1);
r2(DBCIndices) = R0(DBCIndices)+tDBCValues;
r2(FreeIndices) = x;

% Compute gradient
f_r = build_grad_r(atoms,bonds,r2,z);
f_r = f_r(FreeIndices);

% Compute Hessian
[I,J,S] = build_hess_r(atoms,bonds,r2,z);
K_r = sparse(I(:),J(:),S(:),length(r2),length(r2));
K_r = K_r(FreeIndices,FreeIndices);

end
