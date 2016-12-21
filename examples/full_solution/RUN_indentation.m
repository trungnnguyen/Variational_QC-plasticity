%% RUN_indentation: script for example in Section 5.3 (indentation test), full lattice solution

%% Clear workspace
clc; % clear command line
close all; % close all figures
warning('on'); % turn on warnings
matlabrc; % restore MATLAB path, etc.
path([pwd,'/mfiles'],path); % add path to mfiles folder containing all m-files
path([pwd,'/mex'],path); % add path to mex folder containing all mex files

%% Input variables
% Geometry: construct a rectangle [-SizeX,SizeX]x[-SizeY,SizeY]
dSizeX = 1; % lattice spacing along x-axis
dSizeY = 1; % lattice spacing along y-axis
SizeX = 128*dSizeX; % half size of the domain along x-axis
SizeY = 64*dSizeY; % half size of the domain along y-axis

% Material: Potential = [E,H,\sigma_0,\rho]
Potential = [1,10,0.01,0.5];

% Indenter
Rind = SizeY/2; % half-side of the square, or radius of the circular indenter
indenter = 'circle'; % indneter can be 'square' or 'circle'

% Tolerances
TOL_am = 1e-6; % alternating minimization relative tolerance
TOL_r = 1e-6; % elasticity solver relative tolerance
TOL_z = 1e-6; % plasticity solver relative tolerance
TOL_g = 1e-10; % geometric tolerance; a positive number < TOL is treated as zero

%% Construct data structures: atoms and bonds
atoms = build_atoms([SizeX,SizeY,dSizeX,dSizeY]);
R0 = [atoms(:).R]';
bonds = build_bonds(atoms,[SizeX,SizeY,dSizeX,dSizeY,SizeX,SizeY],...
    Potential,Potential,2,TOL_g);

%% Prescribe boundary conditions BC: 0 - no BC, 1 - Dirichlet BC type
% BC = [atomID, 1, 1, uDx, uDy]

% Find atom IDs for individual parts of the boundary
IDGamma_124 = find(abs(R0(1:2:end))>SizeX-TOL_g | R0(2:2:end)<-SizeY+TOL_g); % Gamma_1, Gamma_2, Gamma_4
IDGamma_3 = find(R0(2:2:end)>SizeY-TOL_g & abs(R0(1:2:end))<SizeX-TOL_g); % Gamma_3
IDActive = []; % identifies which atoms from IDGamma_3 are active inequality constraints
lambda = zeros(size(IDGamma_3)); % Lagrange multipliers associated with constraints;

% Add Gamma_124 constraints: no vertical nor horizontal displacements
BC = [IDGamma_124,ones(size(IDGamma_124)),ones(size(IDGamma_124)),...
    zeros(size(IDGamma_124)),zeros(size(IDGamma_124))];

% Extract code numbers
DBCIndices = [2*BC(BC(:,2)==1,1)-1;2*BC(BC(:,3)==1,1)]; % code numbers for constrained atoms
DBCValues = [BC(BC(:,2)==1,4);BC(BC(:,3)==1,5)]; % prescribed displacements for constrained atoms
FreeIndices = setdiff(1:2*length(atoms),DBCIndices)'; % all free code numbers

%% Solve for the evolution process of the system
Time = linspace(0,1,10); % parametrization pseudo-time
fprintf('Solving...\n'),t_start_1 = tic;
fprintf('Time step %d, %d AM it., %d Nwtn it. %g s\n',1,0,0,0);
altIter = zeros(length(Time),1);
nIter = zeros(length(Time),1);
R = [R0,zeros(2*length(atoms),length(Time)-1)]; % r for all time steps
Z = zeros(length(bonds),length(Time)); % z_p for all time steps
K = zeros(length(bonds),length(Time)); % z_c for all time steps
for i = 2:length(Time)
    
    % CENTRE OF THE INDENTER
    C = [0;SizeY+Rind-Time(i)*Rind/4];
    
    % ALTERNATING MINIMIZATION ALGORITHM
    t_start_2 = tic;
    eps_am = TOL_am+1;
    iter = 0;
    Nwtnit = 0;
    r1 = R(:,i-1);
    z1 = Z(:,i-1);
    if norm(z1) == 0; normConstz = 1; else normConstz = norm(z1); end;
    while eps_am > TOL_am;
        
        % MINIMIZE WITH RESPECT TO r
        % [r2,Niter,IDActive,lambda] = minimize_r_I(atoms,bonds,r1,...
        %     DBCIndices,DBCValues,FreeIndices,IDGamma_3,C,Rind,R0,z1,...
        %     TOL_r,TOL_g,indenter,IDActive,lambda);
        
        [r2,Niter,IDActive,lambda] = minimize_r_ISQP(atoms,bonds,r1,...
            DBCIndices,DBCValues,FreeIndices,IDGamma_3,C,Rind,R0,z1,...
            TOL_r,TOL_g,indenter,IDActive,lambda);
        
        % MINIMIZE WITH RESPECT TO z
        z2 = return_mapping(atoms,bonds,r2,Z(:,i-1),K(:,i-1),TOL_z);
        
        % Upadte the AM error
        eps_am = norm(z2-z1)/normConstz;
        
        % Update the results
        r1 = r2;
        z1 = z2;
        iter = iter+1;
        Nwtnit = Nwtnit+Niter;
    end
    
    % Store converged minimizers
    R(:,i) = r2;
    Z(:,i) = z2;
    K(:,i) = K(:,i-1)+abs(Z(:,i)-Z(:,i-1));
    
    % Print iteration message
    fprintf('Time step %d, %d AM it., %d Nwtn it., %g s\n',i,iter,...
        Nwtnit,toc(t_start_2));
    altIter(i) = iter;
    nIter(i) = Nwtnit;
end

% Print AM statistics
fprintf('\nE[AM iter] %g\n',mean(altIter));
fprintf('Max[AM iter] %g\n',max(altIter));
fprintf('E[Nwtn iter] %g\n',mean(nIter));
fprintf('Max[Nwtn iter] %g\n',max(nIter));
fprintf('Time consumed %g s\n',toc(t_start_1));

%% Draw results
% Plot atoms in deformed configuration
scale = 1;
step = length(Time);
figure(1),clf,hold all,axis equal,...
    title(['Deformed structure, scale = ',num2str(scale)]),xlabel('x'),...
    ylabel('y');
plot(R0(1:2:end)+scale*(R(1:2:end,step)-R0(1:2:end)),...
    R0(2:2:end)+scale*(R(2:2:end,step)-R0(2:2:end)),'.k','MarkerSize',5);

% Compute reactions
REACT = zeros(size(R));
for i = 2:length(Time)
    REACT(:,i) = build_grad_r(atoms,bonds,R(:,i),Z(:,i));
end

% Compute energy evolutions
En = zeros(size(Time)); % elastic energy
DissDist = zeros(size(Time)); % dissipation distance
Wext = zeros(1,length(Time)); % work done by the external forces (reactions)
for i = 2:length(Time)
    En(i) = build_en(atoms,bonds,Z(:,i),R(:,i),K(:,i));
    DissDist(i) = build_diss(atoms,bonds,Z(:,i),Z(:,i-1));
    Wext(i) = Wext(i-1)+0.5*(REACT(:,i-1)+REACT(:,i))'*(R(:,i)-R(:,i-1));
end
VarD = cumsum(DissDist); % dissipated energy

% Plot energy evolution paths
figure(2),clf,hold all,box on,xlabel('time [s]'),ylabel('Energy [kJ]'),...
    title('MSPlastInterpolSumPlast');
plot(Time,En+VarD,'--k');
plot(Time,En,'k');
plot(Time,VarD,'-.k');
plot(Time,Wext,':k');
legend('E+Var_D','E','Var_D','W_{ext}','location','northwest');

% Plot vertical reaction force of the indenter
figure(3),clf,hold all,box on,xlabel('time [s]'),...
    ylabel('Reaction');
plot(Time,-sum(REACT(2*IDGamma_3,:)),'k');
legend('vert. react.');

% Deformed state and plastic deformations
id1 = find(abs(R0(1:2:end))<24+TOL_g & R0(2:2:end)>SizeY-24-TOL_g);
id2 = unique([atoms(id1).BondList]);
figure(4); clf, hold all, axis equal,box on,xlabel('x'),ylabel('y');
colormap(jet);
cm = colormap;
minomega = min(K(:,end));
maxomega = max(K(:,end));
for i = 1:length(id2)
    bond = id2(i);
    alpha = bonds(bond).Atoms(1);
    beta = bonds(bond).Atoms(2);
    ra = [R(2*alpha-1,end),R(2*alpha,end)]';
    rb = [R(2*beta-1,end),R(2*beta,end)]';
    omega = K(bond,end);
    colorID = max(1,sum((omega-minomega)/(maxomega-minomega)>...
        0:1/length(cm(:,1)):1));
    colour = cm(colorID,:);
    plot([ra(1),rb(1)],[ra(2),rb(2)],'color',colour);
end
caxis([minomega maxomega]);
colorbar;

%% Call OpenGL postprocessing tool
if exist('draw_OpenGL.exe','file') && exist('freeglut.dll','file') &&...
        exist('glew32.dll','file') && ispc
    call_OpenGl(atoms,[],[],bonds,[],[],R,Z,K,Time,[SizeX,SizeY]);
end
