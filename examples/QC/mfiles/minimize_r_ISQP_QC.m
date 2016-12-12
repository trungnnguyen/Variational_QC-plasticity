function [r2,r,Niter,IDActiveQC,lambda] = minimize_r_ISQP_QC(atoms,...
    samplingatoms,bonds,r1,DBCIndicesQC,tDBCValuesQC,FreeIndicesQC,...
    IDGamma_3QC,Phi,C,Rcut,R0QC,z1,TOL_r,TOL_g,indenter,IDActiveQC,lambda)

% Solve for x, a vector of free degrees of freedom, using primal-dual formulation
r2 = r1;
x = r1(FreeIndicesQC);
Niter = 0;
eps_r = 1+TOL_r;
while eps_r > TOL_r % major iteration
    Niter = Niter+1;
    
    % Taylor expansion of the constraints
    if strcmp(indenter,'circle')
        IA = zeros(2*length(IDGamma_3QC),1);
        JA = zeros(2*length(IDGamma_3QC),1);
        SA = zeros(2*length(IDGamma_3QC),1);
        IH = zeros(2*length(IDGamma_3QC),1);
        SH = zeros(2*length(IDGamma_3QC),1);
        
        % Assembly constraints
        for j = 1:length(IDGamma_3QC)
            ids = [2*IDGamma_3QC(j)-1,...
                2*IDGamma_3QC(j)];
            ra = r2(ids);
            
            % Add gradients of all active inequality constraints
            IA(2*j-1:2*j) = j;
            JA(2*j-1:2*j) = ids;
            SA(2*j-1:2*j) = -2*[ra(1)-C(1),ra(2)-C(2)];
            
            % Build the Hessian of the nonlinear inequality constraints
            IH(2*j-1:2*j) = ids;
            SH(2*j-1:2*j) = -2*lambda(j);
        end
        A = sparse(IA,JA,SA,length(IDGamma_3QC),length(R0QC));
        Hi = sparse(IH,IH,SH,length(R0QC),length(R0QC));
        
        % In the first time step, first iteration, initialize the active set
        if norm(r1-R0QC)==0 && Niter==1
            IDActiveQC = find(Rcut^2-(r2(2*IDGamma_3QC-1)-C(1)).^2-...
                (r2(2*IDGamma_3QC)-C(2)).^2>-TOL_g);
        end
        
        % Violation of the constraints
        B = Rcut^2-(r2(2*IDGamma_3QC-1)-C(1)).^2-...
            (r2(2*IDGamma_3QC)-C(2)).^2;
        
    elseif strcmp(indenter,'square')
        IA = zeros(length(IDGamma_3QC),1);
        JA = zeros(length(IDGamma_3QC),1);
        SA = zeros(length(IDGamma_3QC),1);
        
        % Assembly constraints
        for j = 1:length(IDGamma_3QC)
            
            % Add the condition for the vertical coordinate
            IA(j) = j;
            JA(j) = 2*IDGamma_3QC(j);
            SA(j) = 1;
        end
        A = sparse(IA,JA,SA,length(IDGamma_3QC),length(R0QC));
        Hi = sparse(length(R0QC),length(R0QC));
        
        % In the first time step, first iteration, initialize the active set
        if norm(r1-R0QC)==0 && Niter==1
            IDActiveQC = find(r2(2*IDGamma_3QC)-C(2)+Rcut>-TOL_g &...
                abs(r2(2*IDGamma_3QC-1))<Rcut+TOL_g);
        end
        
        % Violation of the constraints
        B = -Rcut*ones(size(IDGamma_3QC));
        id = find(abs(r2(2*IDGamma_3QC-1))<Rcut+TOL_g);
        B(id) = r2(2*IDGamma_3QC(id))-C(2)+Rcut;
    end
    A = A(:,FreeIndicesQC);
    
    % Taylor expansion of the energy
    [~,G,H] = grad_hess_QC(atoms,samplingatoms,bonds,x,...
        DBCIndicesQC,tDBCValuesQC,FreeIndicesQC,R0QC,z1,Phi);
    H = H + Hi(FreeIndicesQC,FreeIndicesQC);  % add constraint curvatures, based on an estimate of the Lagrange multipliers
    
    % Quadratic programming
    dx1 = zeros(size(x));
    SQPiter = 0;
    while 1 % active set iteration
        SQPiter = SQPiter+1;
        
        % Test optimality conditions
        [minlambda,id0] = min(lambda);
        if minlambda<-TOL_g
            IDActiveQC = setdiff(IDActiveQC,id0);
            lambda(id0) = 0;
        else
            if max(A*dx1+B)<TOL_r && norm(G+H*dx1+A'*lambda)<TOL_r
                break;
            end
        end
        
        % Solve for direction increment
        EH = [H,A(IDActiveQC,:)'
            A(IDActiveQC,:),sparse(length(IDActiveQC),length(IDActiveQC))];
        EG = [G+H*dx1;A(IDActiveQC,:)*dx1+B(IDActiveQC)];
        
        % Solve the system, use primal-dual formulation
        du = -EH\EG;
        dx2 = du(1:end-length(IDActiveQC));
        lambda(IDActiveQC) = du(end-length(IDActiveQC)+1:end); % associated Lagrange multipliers
        
        % Initialize a feasible starting point
        if max(A*dx1+B)>TOL_g && SQPiter==1
            dx1 = dx1+dx2;
            IDActiveQC = find(A*dx1+B>-TOL_g);
            continue;
        end
        
        % Perform line search
        id1 = setdiff((1:length(IDGamma_3QC))',IDActiveQC); % inactive constraints
        id2 = find(A(id1,:)*dx2>0); % those that are about to violate ineq constraints
        if ~isempty(id2)
            [alpha,id3] = min([1;...
                (-A(id1(id2),:)*dx1-B(id1(id2)))./...
                (A(id1(id2),:)*dx2)]); % the one with minimum alpha is taken as active
        else
            alpha = 1;
        end
        dx1 = dx1+alpha*dx2;
        if alpha<1-1e-14
            IDActiveQC = union(IDActiveQC,id1(id2(id3-1)));
        end
        
        % Test for degeneracy
        if norm(dx2)<1e-14
            warning('Degeneracy occured. Stopping SQP.');
            break;
        end
    end
    
    % Update solution
    x = x+dx1;
    
    % Reconstruct converged r2- and r-vectors from x
    r2 = zeros(size(R0QC,1),1);
    r2(DBCIndicesQC) = R0QC(DBCIndicesQC)+tDBCValuesQC;
    r2(FreeIndicesQC) = x; % QC displacement vector reconstruction
    
    % Update the error
    eps_r = norm(dx1)/max(norm(x-R0QC(FreeIndicesQC)),TOL_g)+...
        norm(G+A'*lambda)+norm(max(B,0));
end
r = Phi*r2; % full r reconstruction

end
