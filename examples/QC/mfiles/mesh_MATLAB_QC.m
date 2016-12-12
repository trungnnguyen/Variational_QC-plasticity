function [p,t] = mesh_MATLAB_QC(ExampleType,SizeX,SizeY,dSizeX,dSizeY,...
    FineX,FineY)

%% Input parameters
% Maximal lenght of mesh triangle edge
HMax = max(SizeX,SizeY)/2;

%% Generate input data for initmesh
fprintf('Generating mesh through MATLAB...\n');
tic;

% UNIFORM LOADING TEST
if strcmp(ExampleType,'uniform_loading')
    
    % Create decomposed geometry matrix
    gd = zeros(10,1+FineX/dSizeX*FineY/dSizeY);
    ns = zeros(5,1);
    
    % Upper right rectangle
    X = [0;SizeX;SizeX;0];
    Y = [0;0;SizeY;SizeY];
    SolidData = [3;4;X;Y];
    gd(1:length(SolidData),1) = SolidData;
    sf = 'R1';
    asciiName = double(['R',num2str(1)]);
    ns(1:length(asciiName),1) = asciiName;
    
    % Upper right part of the fully-resolved region
    i = 2;
    for m = 1:FineY/dSizeY
        for n = 1:FineX/dSizeX
            
            % Geometry description matrix
            X = [(n-1)*dSizeX;n*dSizeX;n*dSizeX;(n-1)*dSizeX];
            Y = [(m-1)*dSizeY;(m-1)*dSizeY;m*dSizeY;m*dSizeY];
            SolidData = [3;4;X;Y];
            gd(1:length(SolidData),i) = SolidData;
            
            % Set formula sf
            sf = [sf,'+R',num2str(i)];
            
            % Name space ns
            if i < 10000
                asciiName = double(['R',num2str(i)]);
                ns(1:length(asciiName),i) = asciiName;
            else
                warning('Number of rectangles overflow, max 10,000.');
                return;
            end
            i = i+1;
        end
    end
    
    % PURE BENDING TEST
elseif strcmp(ExampleType,'pure_bending') ||...
        strcmp(ExampleType,'indentation')
    
    % Create decomposed geometry matrix
    gd = zeros(10,1+FineX/dSizeX*FineY/dSizeY);
    ns = zeros(5,1);
    
    % Right rectangle
    X = [0;SizeX;SizeX;0];
    Y = [-SizeY;-SizeY;SizeY;SizeY];
    SolidData = [3;4;X;Y];
    gd(1:length(SolidData),1) = SolidData;
    sf = 'R1';
    asciiName = double(['R',num2str(1)]);
    ns(1:length(asciiName),1) = asciiName;
    
    % Add right half of the fully-resolved region
    i = 2;
    for m = 1:FineY/dSizeY
        for n = 1:FineX/dSizeX
            
            % Geometry description matrix
            X = [(n-1)*dSizeX;n*dSizeX;n*dSizeX;(n-1)*dSizeX];
            Y = [-SizeY+(m-1)*dSizeY;-SizeY+(m-1)*dSizeY;...
                -SizeY+m*dSizeY;-SizeY+m*dSizeY];
            SolidData = [3;4;X;Y];
            gd(1:length(SolidData),i) = SolidData;
            
            % Set formula sf
            sf = [sf,'+R',num2str(i)];
            
            % Name space ns
            if i < 10000
                asciiName = double(['R',num2str(i)]);
                ns(1:length(asciiName),i) = asciiName;
            else
                warning('Number of rectangles overflow, max 10,000.');
                return;
            end
            i = i+1;
        end
    end
end

%% Assembly mesh

% Generate decomposed geometry
dl = decsg(gd,sf,ns);

% Mesh initialization
[p,~,t] = initmesh(dl,'Hmax',HMax,'Hgrad',1.9);
t = t(1:3,:);

% Round coordinates with respect to the underlying lattice
pp = zeros(size(p));
pp(1,:) = round(p(1,:)/dSizeX)*dSizeX;
pp(2,:) = round(p(2,:)/dSizeY)*dSizeY;
tt = t;

% Exclude triangles with zero area
i1 = zeros(size(t,2),1);
k = 0;
for i = 1:size(t,2)
    if polyarea(pp(1,tt(:,i)),pp(2,tt(:,i))) < 1e-3*(dSizeX*dSizeY/2)
        k = k+1;
        i1(k) = i;
    end
end
i1 = i1(1:k);
tt(:,i1) = [];
p = pp;
t = tt;

%% Take mirror images of generated meshes

% UNIFORM LOADING TEST
if strcmp(ExampleType,'uniform_loading')
    
    % Mirror mesh with respect to y-axis
    Mt = t(1:3,:)+size(p,2);
    Mp = [-p(1,:)
        p(2,:)];
    t = [t,Mt];
    p = [p,Mp];
    
    % Mirror the result with respect to x-axis
    Mt = t(1:3,:)+size(p,2);
    Mp = [p(1,:)
        -p(2,:)];
    t = [t,Mt];
    p = [p,Mp];
    
    % PURE BENDING TEST
elseif strcmp(ExampleType,'pure_bending')
    
    % Mirror mesh with respect to y-axis
    Mt = t(1:3,:)+size(p,2);
    Mp = [-p(1,:)
        p(2,:)];
    t = [t,Mt];
    p = [p,Mp];
    
    % INDENTATION
elseif strcmp(ExampleType,'indentation')
    
    % Mirror mesh with respect to y-axis
    Mt = t(1:3,:)+size(p,2);
    Mp = [-p(1,:)
        p(2,:)];
    t = [t,Mt];
    p = [p,Mp];
    
    % Rotate mesh around x-axis
    p = [p(1,:)
        -p(2,:)];
end
fprintf('time consumed %g s\n',toc);

end
