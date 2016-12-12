function [r2,Niter,IDActive,lambdaout] = minimize_r_I(atoms,bonds,r1,...
    DBCIndices,tDBCValues,FreeIndices,IDGamma_3,C,Rcut,R0,z1,TOL_r,...
    TOL_g,indenter,IDActive,lambdain)

% Solve for x, a vector of free degrees of freedom, using primal-dual formulation
x = r1(FreeIndices);
lambda = lambdain(IDActive);
converged = 0;
Niter = 0;
while ~converged
    eps_r = 1+TOL_r;
    while eps_r > TOL_r
        Niter = Niter+1;
        [r2,f_r,K_r] = grad_hess(atoms,bonds,x,DBCIndices,tDBCValues,...
            FreeIndices,R0,z1);
        
        % Introduce indenter through inequality constraints according to IDactive
        if ~isempty(IDActive)
            if strcmp(indenter,'circle')
                IA = zeros(2*length(IDActive),1);
                JA = zeros(2*length(IDActive),1);
                SA = zeros(2*length(IDActive),1);
                IK = zeros(2*length(IDActive),1);
                SK = zeros(2*length(IDActive),1);
                
                % Assembly constraints
                for j = 1:length(IDActive)
                    ids = [2*IDGamma_3(IDActive(j))-1,...
                        2*IDGamma_3(IDActive(j))];
                    ra = r2(ids);
                    
                    % Gradient of active inequality constraints
                    IA(2*j-1:2*j) = j;
                    JA(2*j-1:2*j) = ids;
                    SA(2*j-1:2*j) = -2*[ra(1)-C(1),ra(2)-C(2)];
                    
                    % Hessian of active inequality constraints
                    IK(2*j-1:2*j) = ids;
                    SK(2*j-1:2*j) = -2*lambda(j);
                end
                A = sparse(IA,JA,SA,length(IDActive),length(R0));
                K_i = sparse(IK,IK,SK,length(R0),length(R0));
                K_r = K_r + K_i(FreeIndices,FreeIndices); % add constraint curvatures
                
                % Violation of the constraints
                B = Rcut^2-(r2(2*IDGamma_3(IDActive)-1)-C(1)).^2-...
                    (r2(2*IDGamma_3(IDActive))-C(2)).^2;
            elseif strcmp(indenter,'square')
                IA = zeros(length(IDActive),1);
                JA = zeros(length(IDActive),1);
                SA = zeros(length(IDActive),1);
                
                % Assembly constraints
                for j = 1:length(IDActive)
                    
                    % Add the condition for the vertical coordinate
                    IA(j) = j;
                    JA(j) = 2*IDGamma_3(IDActive(j));
                    SA(j) = 1;
                end
                A = sparse(IA,JA,SA,length(IDActive),length(R0));
                
                % Violation of the constraints
                B = r2(2*IDGamma_3(IDActive))-C(2)+Rcut;
            end
            EK_r = [K_r,A(:,FreeIndices)'
                A(:,FreeIndices),sparse(size(A,1),size(A,1))];
            Ef_r = [f_r;B];
        else
            A = sparse(0,length(R0));
            B = [];
            EK_r = K_r;
            Ef_r = f_r;
        end
        
        % Solve the system, use primal-dual formulation
        du = -EK_r\Ef_r;
        dx = du(1:end-length(IDActive));
        lambda = du(end-length(IDActive)+1:end); % associated Lagrange multipliers
        x = x+dx;
        
        % Update the error
        eps_r = norm(dx)/max(norm(x-R0(FreeIndices)),TOL_g)+...
            norm(f_r+A(:,FreeIndices)'*lambda)+norm(B);
        if norm(dx)<TOL_g && norm(x-R0(FreeIndices))<TOL_g
            eps_r = 0;
            warning('Division by 0, taking eps_r = 0');
        end
    end
    
    % Reconstruct the converged r-vector from the iterative variable x
    r2(DBCIndices) = R0(DBCIndices)+tDBCValues;
    r2(FreeIndices) = x;
    
    % Test for inequality constraints, update IDActive
    if strcmp(indenter,'circle')
        tlambda = zeros(size(IDGamma_3));
        tlambda(IDActive) = lambda;
        tIDActive = find(Rcut^2-(r2(2*IDGamma_3-1)-C(1)).^2-...
            (r2(2*IDGamma_3)-C(2)).^2>-TOL_g);
        
        % Test Lagrange multipliers
        tIDActive = setdiff(tIDActive,IDActive(lambda<0));
        lambda = tlambda(tIDActive);
    elseif strcmp(indenter,'square')
        tIDActive = find(r2(2*IDGamma_3)-C(2)+Rcut>-TOL_g &...
            abs(r2(2*IDGamma_3-1))<Rcut+TOL_g);
        
        % Test Lagrange multipliers
        tIDActive = setdiff(tIDActive,IDActive(lambda<0));
    end
    if isempty(union(setdiff(tIDActive,IDActive),...
            setdiff(IDActive,tIDActive)))
        converged = 1;
    end
    IDActive = tIDActive;
end
lambdaout = zeros(size(IDGamma_3));
lambdaout(IDActive) = lambda;

end
