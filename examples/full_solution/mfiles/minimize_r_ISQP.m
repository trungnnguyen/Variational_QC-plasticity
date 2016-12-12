function [r2,Niter,IDActive,lambda] = minimize_r_ISQP(atoms,bonds,r1,...
    DBCIndices,tDBCValues,FreeIndices,IDGamma_3,C,Rcut,R0,z1,TOL_r,...
    TOL_g,indenter,IDActive,lambda)

% Solve for x, a vector of free degrees of freedom, using primal-dual formulation
r2 = r1;
x = r1(FreeIndices);
Niter = 0;
eps_r = 1+TOL_r;
while eps_r > TOL_r % major iteration
    Niter = Niter+1;
    
    % Taylor expansion of the constraints
    if strcmp(indenter,'circle')
        IA = zeros(2*length(IDGamma_3),1);
        JA = zeros(2*length(IDGamma_3),1);
        SA = zeros(2*length(IDGamma_3),1);
        IK = zeros(2*length(IDGamma_3),1);
        SK = zeros(2*length(IDGamma_3),1);
        
        % Assembly constraints
        for j = 1:length(IDGamma_3)
            ids = [2*IDGamma_3(j)-1,...
                2*IDGamma_3(j)];
            ra = r2(ids);
            
            % Add gradients of all active inequality constraints
            IA(2*j-1:2*j) = j;
            JA(2*j-1:2*j) = ids;
            SA(2*j-1:2*j) = -2*[ra(1)-C(1),ra(2)-C(2)];
            
            % Build the Hessian of the nonlinear inequality constraints
            IK(2*j-1:2*j) = ids;
            SK(2*j-1:2*j) = -2*lambda(j);
        end
        A = sparse(IA,JA,SA,length(IDGamma_3),length(R0));
        Ki = sparse(IK,IK,SK,length(R0),length(R0));
        
        % In the first time step, first iteration, initialize the active set
        if norm(r1-R0)==0 && Niter==1
            IDActive = find(Rcut^2-(r2(2*IDGamma_3-1)-C(1)).^2-...
                (r2(2*IDGamma_3)-C(2)).^2>-TOL_g);
        end
        
        % Violation of the constraints
        B = Rcut^2-(r2(2*IDGamma_3-1)-C(1)).^2-...
            (r2(2*IDGamma_3)-C(2)).^2;
        
    elseif strcmp(indenter,'square')
        IA = zeros(length(IDGamma_3),1);
        JA = zeros(length(IDGamma_3),1);
        SA = zeros(length(IDGamma_3),1);
        
        % Assembly constraints
        for j = 1:length(IDGamma_3)
            
            % Add the condition for the vertical coordinate
            IA(j) = j;
            JA(j) = 2*IDGamma_3(j);
            SA(j) = 1;
        end
        A = sparse(IA,JA,SA,length(IDGamma_3),length(R0));
        Ki = sparse(length(R0),length(R0));
        
        % In the first time step, first iteration, initialize the active set
        if norm(r1-R0)==0 && Niter==1
            IDActive = find(r2(2*IDGamma_3)-C(2)+Rcut>-TOL_g &...
                abs(r2(2*IDGamma_3-1))<Rcut+TOL_g);
        end
        
        % Violation of the constraints
        B = -Rcut*ones(size(IDGamma_3));
        id = find(abs(r2(2*IDGamma_3-1))<Rcut+TOL_g);
        B(id) = r2(2*IDGamma_3(id))-C(2)+Rcut;
    end
    A = A(:,FreeIndices);
    
    % Taylor expansion of the energy
    [~,f_r,K_r] = grad_hess(atoms,bonds,x,DBCIndices,tDBCValues,...
        FreeIndices,R0,z1);
    K_r = K_r + Ki(FreeIndices,FreeIndices); % add constraint curvatures, based on an estimate of the Lagrange multipliers
    
    % Quadratic programming
    dx1 = zeros(size(x));
    SQPiter = 0;
    while 1 % active set iteration
        SQPiter = SQPiter+1;
        
        % Test optimality conditions
        [minlambda,id0] = min(lambda);
        if minlambda<-TOL_g
            IDActive = setdiff(IDActive,id0);
            lambda(id0) = 0;
        else
            if max(A*dx1+B)<TOL_r && norm(f_r+K_r*dx1+A'*lambda)<TOL_r
                break;
            end
        end
        
        % Solve for direction increment
        EK = [K_r,A(IDActive,:)'
            A(IDActive,:),sparse(length(IDActive),length(IDActive))];
        Ef = [f_r+K_r*dx1;A(IDActive,:)*dx1+B(IDActive)];
        
        % Solve the system, use primal-dual formulation
        du = -EK\Ef;
        dx2 = du(1:end-length(IDActive));
        lambda(IDActive) = du(end-length(IDActive)+1:end); % associated Lagrange multipliers
        
        % Initialize a feasible starting point
        if max(A*dx1+B)>TOL_g && SQPiter==1
            dx1 = dx1+dx2;
            IDActive = find(A*dx1+B>-TOL_g);
            continue;
        end
        
        % Perform line search
        id1 = setdiff((1:length(IDGamma_3))',IDActive); % inactive constraints
        id2 = find(A(id1,:)*dx2>0); % those that are about to violate ineq constraints
        if ~isempty(id2)
            [alpha,id3] = min([1;...
                (-A(id1(id2),:)*dx1-B(id1(id2)))./...
                (A(id1(id2),:)*dx2)]); % the one with minimum alpha is taken as active
        else
            alpha = 1;
        end
        dx1 = dx1+alpha*dx2;
        if alpha<1-1e-14 % && alpha~=0
            IDActive = union(IDActive,id1(id2(id3-1)));
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
    r2 = zeros(size(R0,1),1);
    r2(DBCIndices) = R0(DBCIndices)+tDBCValues;
    r2(FreeIndices) = x;
    
    % Update the error
    eps_r = norm(dx1)/max(norm(x-R0(FreeIndices)),TOL_g)+...
        norm(f_r+A'*lambda)+norm(max(B,0));
end

end
