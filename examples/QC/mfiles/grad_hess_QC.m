function [r2,G,H] = grad_hess_QC(atoms,samplingatoms,bonds,x,...
    DBCIndicesQC,tDBCValuesQC,FreeIndicesQC,R0QC,z,Phi)

% Reconstruct r
r2 = zeros(size(R0QC));
r2(DBCIndicesQC) = R0QC(DBCIndicesQC)+tDBCValuesQC;
r2(FreeIndicesQC) = x;
r = Phi*r2;

% Compute gradient
f_r = build_grad_r_QC(atoms,samplingatoms,bonds,r,z);
G = Phi'*f_r;
G = G(FreeIndicesQC);

% Compute Hessian
[I,J,S] = build_hess_r_QC(atoms,samplingatoms,bonds,r,z);
K_r = sparse(I(:),J(:),S(:),length(r),length(r));
H = Phi'*K_r*Phi;
H = H(FreeIndicesQC,FreeIndicesQC);

end
