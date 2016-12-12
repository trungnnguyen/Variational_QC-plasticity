%% RUN_pure_bending: script for example in Section 5.2 (pure bending test), full lattice

%% Clear workspace
clc; % clear command line
close all; % close all figures
matlabrc; % restore MATLAB path, etc.
path([pwd,'/mfiles'],path); % add path to mfiles folder containing all m-files
path([pwd,'/mex'],path); % add path to mex folder containing all mex files

%% Input variables
% Geometry: construct a rectangle [-SizeX,SizeX]x[-SizeY,SizeY]
dSizeX = 1; % lattice spacing along x-axis
dSizeY = 1; % lattice spacing along y-axis
SizeX = 100*dSizeX; % half size of the domain along x-axis
SizeY = 50*dSizeY; % half size of the domain along y-axis
RigidX = 3*dSizeX; % half size of the bottom-edge inclusion along x-axis
RigidY = 6*dSizeY; % size of the bottom-edge inclusion along y-axis

% Material: Potential = [E,H,\sigma_0,\rho]
PotentialI = [1*100,10*100,0.01*1000,0.5]; % for Inclusion
PotentialM = [1,10,0.01,0.5]; % for surrounding Matrix

% Tolerances
TOL_am = 1e-6; % alternating minimization relative tolerance
TOL_r = 1e-6; % elasticity solver relative tolerance
TOL_z = 1e-6; % plasticity solver relative tolerance
TOL_g = 1e-10; % geometric tolerance; a positive number < TOL is treated as zero

%% Construct data structures: atoms and bonds
atoms = build_atoms([SizeX,SizeY,dSizeX,dSizeY]);
R0 = [atoms(:).R]';
bonds = build_bonds(atoms,[SizeX,SizeY,dSizeX,dSizeY,RigidX,RigidY],...
    PotentialI,PotentialM,1,TOL_g);

%% Prescribe boundary conditions BC: 0 - no BC, 1 - Dirichlet BC type
% BC = [atomID, 1, 1, uDx, uDy]
THETA = pi/6; % target bending angle

% Find atom IDs for individual parts of the boundary
IDGamma_2 = find(R0(1:2:end)>SizeX-TOL_g); % Gamma_2
IDGamma_4 = find(R0(1:2:end)<-SizeX+TOL_g); % Gamma_4

% Add Gamma_4 constraints: no horizontal displacement - symmetry conditions
BC = [IDGamma_4,ones(size(IDGamma_4)),zeros(size(IDGamma_4)),...
    zeros(size(IDGamma_4)),zeros(size(IDGamma_4))];

% Add fixed corner node
BC(1,3) = 1; % left bottom corner

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
    
    % UPDATE CONSTRAINTS
    theta = THETA*Time(i); % angle theta
    
    % Assembly constraints
    C = zeros(size(IDGamma_2,1)-1,2*length(atoms));
    for j = 1:size(IDGamma_2,1)-1
        
        % Add the slope condition
        C(j,2*IDGamma_2(j+1)-1) = 1;
        C(j,2*IDGamma_2(j+1)) = -tan(theta);
        C(j,2*IDGamma_2(j)-1) = -1;
        C(j,2*IDGamma_2(j)) = tan(theta);
    end
    C = sparse(C);
    
    % ALTERNATING MINIMIZATION
    t_start_2 = tic;
    eps_am = TOL_am+1;
    iter = 0;
    Nwtnit = 0;
    r1 = R(:,i-1);
    z1 = Z(:,i-1);
    if norm(z1) == 0; normConstz = 1; else normConstz = norm(z1); end;
    while eps_am > TOL_am;
        
        % MINIMIZE WITH RESPECT TO r
        [r2,Niter] = minimize_r(atoms,bonds,r1,DBCIndices,...
            Time(i)*DBCValues,FreeIndices,R0,z1,C,TOL_r);
        
        % MINIMIZE WITH RESPECT TO z
        z2 = return_mapping(atoms,bonds,r2,Z(:,i-1),K(:,i-1),TOL_z);
        
        % Upadte the AM error
        eps_am = norm(z2-z1)/normConstz;
        
        % Update results
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
DissDist = zeros(size(Time)); % dissipation distances
Wext = zeros(1,length(Time)); % work done by the external forces (reactions)
for i = 2:length(Time)
    En(i) = build_en(atoms,bonds,Z(:,i),R(:,i),K(:,i));
    DissDist(i) = build_diss(atoms,bonds,Z(:,i),Z(:,i-1));
    Wext(i) = Wext(i-1)+0.5*(REACT(:,i-1)+REACT(:,i))'*(R(:,i)-R(:,i-1));
end
VarD = cumsum(DissDist); % dissipated energy

% Plot energy evolution paths
figure(2),clf,hold all,box on,xlabel('time [s]'),ylabel('Energy [kJ]'),...
    title('Energy evolutions');
plot(Time,En+VarD,'--k');
plot(Time,En,'k');
plot(Time,VarD,'-.k');
plot(Time,Wext,':k');
legend('E+Var_D','E','Var_D','W_{ext}','location','northwest');

%% Call OpenGL postprocessing tool
if exist('draw_OpenGL.exe','file') && exist('freeglut.dll','file') &&...
        exist('glew32.dll','file') && ispc
    call_OpenGl(atoms,[],[],bonds,[],[],R,Z,K,Time,[SizeX,SizeY]);
end
