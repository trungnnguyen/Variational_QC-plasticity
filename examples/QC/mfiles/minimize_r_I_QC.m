function [r2,r,Niter,IDActiveQC,lambdaout] = minimize_r_I_QC(atoms,...
    samplingatoms,bonds,r1,DBCIndicesQC,tDBCValuesQC,FreeIndicesQC,...
    IDGamma_3QC,Phi,C,Rcut,R0QC,z1,TOL_r,TOL_g,indenter,IDActiveQC,lambdain)

% Solve for x, a vector of free degrees of freedom, using primal-dual formulation
x = r1(FreeIndicesQC);
lambda = lambdain(IDActiveQC);
converged = 0;
Niter = 0;
while ~converged
    eps_r = 1+TOL_r;
    while eps_r > TOL_r
        Niter = Niter+1;
        [r2,G,H] = grad_hess_QC(atoms,samplingatoms,bonds,x,...
            DBCIndicesQC,tDBCValuesQC,FreeIndicesQC,R0QC,z1,Phi);
        
        % Introduce indenter through inequality constraints according to IDactive
        if ~isempty(IDActiveQC)
            if strcmp(indenter,'circle')
                IA = zeros(2*length(IDActiveQC),1);
                JA = zeros(2*length(IDActiveQC),1);
                SA = zeros(2*length(IDActiveQC),1);
                IH = zeros(2*length(IDActiveQC),1);
                SH = zeros(2*length(IDActiveQC),1);
                
                % Assembly constraints
                for j = 1:length(IDActiveQC)
                    ids = [2*IDGamma_3QC(IDActiveQC(j))-1,...
                        2*IDGamma_3QC(IDActiveQC(j))];
                    ra = r2(ids);
                    
                    % Add gradients of all active inequality constraints
                    IA(2*j-1:2*j) = j;
                    JA(2*j-1:2*j) = ids;
                    SA(2*j-1:2*j) = -2*[ra(1)-C(1),ra(2)-C(2)];
                    
                    % Build the Hessian of the nonlinear inequality constraints
                    IH(2*j-1:2*j) = ids;
                    SH(2*j-1:2*j) = -2*lambda(j);
                end
                A = sparse(IA,JA,SA,length(IDActiveQC),length(R0QC));
                Hi = sparse(IH,IH,SH,length(R0QC),length(R0QC));
                H = H + Hi(FreeIndicesQC,FreeIndicesQC); % add constraint curvatures
                
                % Violation of the constraints
                B = Rcut^2-(r2(2*IDGamma_3QC(IDActiveQC)-1)-C(1)).^2-...
                    (r2(2*IDGamma_3QC(IDActiveQC))-C(2)).^2;
            elseif strcmp(indenter,'square')
                IA = zeros(length(IDActiveQC),1);
                JA = zeros(length(IDActiveQC),1);
                SA = zeros(length(IDActiveQC),1);
                
                % Assembly constraints
                for j = 1:length(IDActiveQC)
                    
                    % Add the condition for the vertical coordinate
                    IA(j) = j;
                    JA(j) = 2*IDGamma_3QC(IDActiveQC(j));
                    SA(j) = 1;
                end
                A = sparse(IA,JA,SA,length(IDActiveQC),length(R0QC));
                
                % Violation of the constraints
                B = r2(2*IDGamma_3QC(IDActiveQC))-C(2)+Rcut;
            end
            EH = [H,A(:,FreeIndicesQC)'
                A(:,FreeIndicesQC),sparse(size(A,1),size(A,1))];
            EG = [G;B];
        else
            A = sparse(0,length(R0QC));
            B = [];
            EH = H;
            EG = G;
        end
        
        % Solve the system, use primal-dual formulation
        du = -EH\EG;
        dx = du(1:end-length(IDActiveQC));
        lambda = du(end-length(IDActiveQC)+1:end); % associated Lagrange multipliers
        x = x+dx;
        
        % Update the error
        eps_r = norm(dx)/max(norm(x-R0QC(FreeIndicesQC)),TOL_g)+...
            norm(G+A(:,FreeIndicesQC)'*lambda)+norm(B);
        if norm(dx)<TOL_g && norm(x-R0QC(FreeIndicesQC))<TOL_g
            eps_r = 0;
            warning('Division by 0, taking eps_r = 0');
        end
    end
    
    % Reconstruct converged r2- and r-vectors from x
    r2 = zeros(size(R0QC,1),1);
    r2(DBCIndicesQC) = R0QC(DBCIndicesQC)+tDBCValuesQC;
    r2(FreeIndicesQC) = x; % QC displacement vector reconstruction
    
    % Test for inequality constraints, update IDActive
    if strcmp(indenter,'circle')
        tlambda = zeros(size(IDGamma_3QC));
        tlambda(IDActiveQC) = lambda;
        tIDActiveQC = find(Rcut^2-(r2(2*IDGamma_3QC-1)-C(1)).^2-...
            (r2(2*IDGamma_3QC)-C(2)).^2>-TOL_g);
        
        % Test Lagrange multipliers
        tIDActiveQC = setdiff(tIDActiveQC,IDActiveQC(lambda<0));
        lambda = tlambda(tIDActiveQC);
    elseif strcmp(indenter,'square')
        tIDActiveQC = find(r2(2*IDGamma_3QC)-C(2)+Rcut>-TOL_g &...
            abs(r2(2*IDGamma_3QC-1))<Rcut+TOL_g);
        
        % Test Lagrange multipliers
        tIDActiveQC = setdiff(tIDActiveQC,IDActiveQC(lambda<0));
    end
    if isempty(union(setdiff(tIDActiveQC,IDActiveQC),...
            setdiff(IDActiveQC,tIDActiveQC)))
        converged = 1;
    end
    IDActiveQC = tIDActiveQC;
end
r = Phi*r2; % full r reconstruction
lambdaout = zeros(size(IDGamma_3QC));
lambdaout(IDActiveQC) = lambda;

end
