%% RUN_pure_bending_QC: script for example in Section 5.2 (pure bending test), QC solution

%% Clear workspace
clc; % clear command line
close all; % close all figures
warning('on'); % turn on warnings
matlabrc; % restore MATLAB path, etc.
path([pwd,'/mfiles'],path); % add path to mfiles folder containing all m-files
path([pwd,'/mex'],path); % add path to mex folder containing all mex files
path([pwd '/mesh'],path); % add path to folder containing meshes

%% Input variables
% Geometry: construct a rectangle [-SizeX,SizeX]x[-SizeY,SizeY]
dSizeX = 1; % lattice spacing along x-axis
dSizeY = 1; % lattice spacing along y-axis
SizeX = 100*dSizeX; % half size of the domain along x-axis
SizeY = 50*dSizeY; % half size of the domain along y-axis
RigidX = 3*dSizeX; % half size of the central inclusion along x-axis
RigidY = 6*dSizeY; % size of the central inclusion along y-axis

% Material: Potential = [E,H,\sigma_0,\rho]
PotentialI = [1*100,10,0.01*1000,0.5]; % for Inclusion
PotentialM = [1,10,0.01,0.5]; % for surrounding Matrix

% Tolerances
TOL_am = 1e-6; % alternating minimization relative tolerance
TOL_r = 1e-6; % elasticity solver relative tolerance
TOL_z = 1e-6; % plasticity solver relative tolerance
TOL_g = 1e-10; % geometric tolerance; a positive number < TOL is treated as zero

% Summation rule: 0 - exact Beex, 1 - central summation Beex, 2 - full summation
SumRule = 1;

% Select a mesh: 0 - load mesh, 1 - RTIN mesh, 2 - construct T3D mesh, 3 - construct MATLAB mesh
MeshType = 0;

%% Construct data structures: atoms, bonds
atoms = build_atoms([SizeX,SizeY,dSizeX,dSizeY]);
R0 = [atoms(:).R]';
bonds = build_bonds(atoms,[SizeX,SizeY,dSizeX,dSizeY,RigidX,RigidY],...
    PotentialI,PotentialM,1,TOL_g);

%% Interpolation step
FineX = 7*dSizeX; % size of a half of the region along x-axis that will be fully refined
FineY = 10*dSizeY; % size of the region along y-axis that will be fully refined
ExampleType = 'pure_bending';
switch MeshType
    case 0
        fprintf('Load mesh...\n');
        load('Ba','p','t');
    case 1
        [p,t] = mesh_RTIN_QC(ExampleType,SizeX,SizeY,dSizeX,dSizeY,...
            FineX,FineY,TOL_g);
    case 2
        if exist('T3d.exe','file')
            [p,t] = mesh_T3d_QC(ExampleType,SizeX,SizeY,dSizeX,dSizeY,...
                FineX,FineY);
        end
    case 3
        [p,t] = mesh_MATLAB_QC(ExampleType,SizeX,SizeY,dSizeX,dSizeY,...
            FineX,FineY);
end

% Construct triangles, repatoms, R0QC, and Phi for QC
fprintf('Building triangles database and Phi matrix...\n');
tic;
[triangles,repatoms,I,J,S] = sort_atoms_QC(p,t,atoms,TOL_g);
Phi = accumarray([I,J],S,[2*length(atoms),2*length(repatoms)],@max,[],true);
R0QC = [atoms(repatoms).R]';

% Check for hanging nodes
for i = 1:length(triangles)
    if length(triangles(i).NeighTriangles)<3
        % Test that given triangle is boundary one
        c1 = (triangles(i).P1+triangles(i).P2)/2;
        c2 = (triangles(i).P2+triangles(i).P3)/2;
        c3 = (triangles(i).P3+triangles(i).P1)/2;
        if (abs(c1(1))<SizeX-TOL_g && abs(c1(2))<SizeY-TOL_g) &&...
                (abs(c2(1))<SizeX-TOL_g && abs(c2(2))<SizeY-TOL_g) &&...
                (abs(c3(1))<SizeX-TOL_g && abs(c3(2))<SizeY-TOL_g)
            error('Mesh error: hanging nodes detected.');
        end
    end
end
fprintf('time consumed %g s\n',toc);

%% Determine sampling atoms within triangles
fprintf('Building sampling atoms database...\n'),tic;
samplingatoms = sort_sampling_atoms_QC(atoms,triangles,SumRule,...
    [SizeX,SizeY],1,TOL_g);
samplingbondsID = unique([atoms([samplingatoms(:).ID]).BondList])';
fprintf('time consumed %g s\n',toc);

% Plot atoms, repatoms, sampling atoms, and triangulation
figure(1); clf, hold all, axis equal, xlabel('x'), ylabel('y');
pdemesh(p,[],t);
plot(R0(1:2:end),R0(2:2:end),'.k','MarkerSize',1);
plot(R0(2*[samplingatoms(:).ID]-1),R0(2*[samplingatoms(:).ID]),...
    '.m','MarkerSize',5);
plot(R0QC(1:2:end),R0QC(2:2:end),'.g','MarkerSize',5);

%% Prescribe boundary conditions BCQC: 0 - no BC, 1 - Dirichlet BC type
% BCQC = [repatomID, 1, 1, uDx, uDy]
THETA = pi/6; % target bending angle

% Find repatom IDs for individual parts of the boundary
IDGamma_2QC = find(R0QC(1:2:end)>SizeX-TOL_g); % Gamma_2
IDGamma_4QC = find(R0QC(1:2:end)<-SizeX+TOL_g); % Gamma_4

% Add Gamma_4 constraints: no horizontal displacement - symmetry conditions
BCQC = [IDGamma_4QC,ones(size(IDGamma_4QC)),zeros(size(IDGamma_4QC)),...
    zeros(size(IDGamma_4QC)),zeros(size(IDGamma_4QC))];

% Add fixed corner node
BCQC(1,3) = 1; % left bottom corner

% Extract code numbers
DBCIndicesQC = [2*BCQC(BCQC(:,2)==1,1)-1;2*BCQC(BCQC(:,3)==1,1)]; % code numbers for constrained repatoms
DBCValuesQC = [BCQC(BCQC(:,2)==1,4);BCQC(BCQC(:,3)==1,5)]; % prescribed displacements for constrained repatoms
FreeIndicesQC = setdiff(1:2*length(repatoms),DBCIndicesQC)'; % all free code numbers

%% Solve for the evolution process of the system
Time = linspace(0,1,10); % parametrization pseudo-time
fprintf('\nSolving...\n'),t_start_1 = tic;
fprintf('Time step %d, %d AM it., %d Nwtn it. %g s\n',1,0,0,0);
altIter = zeros(length(Time),1);
nIter = zeros(length(Time),1);
R = [R0,zeros(2*length(atoms),length(Time)-1)]; % r for all time steps
RQC = [R0QC,zeros(2*length(repatoms),length(Time)-1)]; % r_rep for all time steps
Z = sparse(length(bonds),length(Time)); % z_p for all time steps
K = sparse(length(bonds),length(Time)); % z_c for all time steps
for i = 2:length(Time)
    
    % UPDATE CONSTRAINTS
    theta = THETA*Time(i); % angle theta
    
    % Assembly constraints
    C = zeros(size(IDGamma_2QC,1)-1,2*length(repatoms));
    for j = 1:size(IDGamma_2QC,1)-1
        
        % Add the slope condition
        C(j,2*IDGamma_2QC(j+1)-1) = 1;
        C(j,2*IDGamma_2QC(j+1)) = -tan(theta);
        C(j,2*IDGamma_2QC(j)-1) = -1;
        C(j,2*IDGamma_2QC(j)) = tan(theta);
    end
    C = sparse(C);
    
    % ALTERNATING MINIMIZATION ALGORITHM
    t_start_2 = tic;
    eps_am = TOL_am+1;
    iter = 0;
    Nwtnit = 0;
    r1 = RQC(:,i-1);
    z1 = full(Z(:,i-1));
    if norm(z1) == 0; normConstz = 1; else normConstz = norm(z1); end;
    while eps_am > TOL_am;
        
        % MINIMIZE WITH RESPECT TO r
        [r2,r,Niter] = minimize_r_QC(atoms,samplingatoms,bonds,r1,...
            DBCIndicesQC,Time(i)*DBCValuesQC,FreeIndicesQC,R0QC,z1,C,...
            Phi,TOL_r);
        
        % MINIMIZE WITH RESPECT TO z
        z2 = return_mapping_QC(atoms,bonds,samplingbondsID,r,...
            full(Z(:,i-1)),full(K(:,i-1)),TOL_z);
        
        % Upadte AM error
        eps_am = norm(z2-z1)/normConstz;
        
        % Update results
        r1 = r2;
        z1 = z2;
        iter = iter+1;
        Nwtnit = Nwtnit+Niter;
    end
    
    % Store converged minimizers
    R(:,i) = r;
    RQC(:,i) = r2;
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
figure(2),clf,hold all,axis equal,title(['Deformed structure, scale = ',...
    num2str(scale)]),xlabel('x'),ylabel('y');
plot(R0(1:2:end)+scale*(R(1:2:end,step)-R0(1:2:end)),...
    R0(2:2:end)+scale*(R(2:2:end,step)-R0(2:2:end)),'.k','MarkerSize',5);

% Compute reactions
REACT = zeros(size(R));
for i = 2:length(Time)
    REACT(:,i) = build_grad_r_QC(atoms,samplingatoms,bonds,...
        R(:,i),full(Z(:,i)));
end
REACTQC = Phi'*REACT;

% Compute energy evolutions
En = zeros(size(Time)); % elastic energy
DissDist = zeros(size(Time)); % dissipation distances
Wext = zeros(1,length(Time)); % work done by external forces (reactions)
for i = 2:length(Time)
    En(i) = build_en_QC(atoms,samplingatoms,bonds,full(Z(:,i)),...
        R(:,i),full(K(:,i)));
    DissDist(i) = build_diss_QC(atoms,samplingatoms,bonds,full(Z(:,i)),...
        full(Z(:,i-1)));
    Wext(i) = Wext(i-1)+0.5*(REACTQC(:,i-1)+...
        REACTQC(:,i))'*(RQC(:,i)-RQC(:,i-1));
end
VarD = cumsum(DissDist); % dissipated energy

% Plot energy evolution paths
figure(3),clf,hold all,box on,xlabel('time [s]'),...
    ylabel('Energy [kJ]'),title('MSPlastInterpolSumPlast');
plot(Time,En+VarD,'--k');
plot(Time,En,'k');
plot(Time,VarD,'-.k');
plot(Time,Wext,':k');
legend('E+Var_D','E','Var_D','W_{ext}','location','northwest');

%% Call OpenGL postprocessing tool
if exist('draw_OpenGL.exe','file') && exist('freeglut.dll','file') &&...
        exist('glew32.dll','file') && ispc
    call_OpenGl(atoms,samplingatoms,repatoms,bonds,samplingbondsID,...
        triangles,R,full(Z),full(K),Time,[SizeX,SizeY]);
end
