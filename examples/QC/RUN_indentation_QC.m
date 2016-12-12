%% RUN_indentation_QC: script for example in Section 5.3 (indentation test), QC solution

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
SizeX = 128*dSizeX; % half size of the domain along x-axis
SizeY = 64*dSizeY; % half size of the domain along y-axis

% Material: Potential = [E,H,\sigma_0,\rho]
Potential = [1,10,0.01,0.5];

% Indenter
Rcut = SizeY/2; % half-side of the square and radius of the circular indenter
indenter = 'circle'; % indneter can be 'square' or 'circle'

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
bonds = build_bonds(atoms,[SizeX,SizeY,dSizeX,dSizeY,SizeX,SizeY],...
    Potential,Potential,2,TOL_g);

%% Interpolation step
FineX = ceil(0.5*Rcut/dSizeX)*dSizeX; % size of a half of the region along x-axis that will be fully refined
FineY = ceil(0.4*Rcut/dSizeY)*dSizeY; % size of the region along y-axis that will be fully refined
ExampleType = 'indentation';
switch MeshType
    case 0
        fprintf('Load mesh...\n');
        load('Ia','p','t');
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
        % Triangle must be a boundary one, otherwise a hanging node exists
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
    [SizeX,SizeY],0,TOL_g);
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

% Find atom IDs for individual parts of the boundary
IDGamma_124QC = find(abs(R0QC(1:2:end))>SizeX-TOL_g | R0QC(2:2:end)<-SizeY+TOL_g); % Gamma_1, Gamma_2, Gamma_4
IDGamma_3QC = find(R0QC(2:2:end)>SizeY-TOL_g & abs(R0QC(1:2:end))<SizeX-TOL_g); % Gamma_3
IDActiveQC = []; % identifies which repatoms from IDGamma_3QC are active inequality constraints
lambda = zeros(size(IDGamma_3QC)); % lagrange multipliers associated with active constraints

% Add Gamma_124 constraints: no vertical nor horizontal displacements
BCQC = [IDGamma_124QC,ones(size(IDGamma_124QC)),ones(size(IDGamma_124QC)),...
    zeros(size(IDGamma_124QC)),zeros(size(IDGamma_124QC))];

% Extract code numbers
DBCIndicesQC = [2*BCQC(BCQC(:,2)==1,1)-1;2*BCQC(BCQC(:,3)==1,1)]; % code numbers for constrained atoms
DBCValuesQC = [BCQC(BCQC(:,2)==1,4);BCQC(BCQC(:,3)==1,5)]; % prescribed displacements for constrained atoms
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
    
    % CENTRE OF THE INDENTER
    C = [0;SizeY+Rcut-Time(i)*Rcut/4];
    
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
        % [r2,r,Niter,IDActiveQC,lambda] = minimize_r_I_QC(atoms,...
        %     samplingatoms,bonds,r1,DBCIndicesQC,DBCValuesQC,...
        %     FreeIndicesQC,IDGamma_3QC,Phi,C,Rcut,R0QC,z1,TOL_r,TOL_g,...
        %     indenter,IDActiveQC,lambda);
        
        [r2,r,Niter,IDActiveQC,lambda] = minimize_r_ISQP_QC(atoms,...
            samplingatoms,bonds,r1,DBCIndicesQC,DBCValuesQC,...
            FreeIndicesQC,IDGamma_3QC,Phi,C,Rcut,R0QC,z1,TOL_r,TOL_g,...
            indenter,IDActiveQC,lambda);
        
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
figure(2),clf,hold all,axis equal,...
    title(['Deformed structure, scale = ',num2str(scale)]),...
    xlabel('x'),ylabel('y');
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
    DissDist(i) = build_diss_QC(atoms,samplingatoms,bonds,...
        full(Z(:,i)),full(Z(:,i-1)));
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

% Plot vertical reaction force of the indenter
figure(4),clf,hold all,box on,xlabel('time [s]'),...
    ylabel('Reaction');
plot(Time,-sum(REACTQC(2*IDGamma_3QC,:)),'k');
legend('vert. react.');

% Deformed state and plastic deformations
id1 = find(abs(R0QC(1:2:end))<24+TOL_g & R0QC(2:2:end)>SizeY-24-TOL_g);
id2 = unique([atoms(repatoms(id1)).BondList]);
figure(4); clf, hold all, axis equal,box on,xlabel('x'),ylabel('y');
colormap(jet);
cm = colormap;
minomega = full(min(K(:,end)));
maxomega = full(max(K(:,end)));
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
    call_OpenGl(atoms,samplingatoms,repatoms,bonds,samplingbondsID,...
        triangles,R,full(Z),full(K),Time,[SizeX,SizeY]);
end
