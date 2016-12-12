function [p,t] = mesh_RTIN_QC(ExampleType,SizeX,SizeY,dSizeX,dSizeY,...
    FineX,FineY,TOL_g)

%% Input parameters
if dSizeX ~= dSizeY
    error('mesh_RTIN_QC: dSizeX must equal dSizeY.');
end

% Maximal lenght of mesh triangle edge
HMax = max(SizeX,SizeY);

%% Build right half of the regular mesh
fprintf('Generating RTIN mesh...\n');
tic;

edgeSize2 = dSizeX;
failure = 1;
while failure
    if mod(SizeX,edgeSize2) == 0 && mod(2*SizeY,edgeSize2) == 0
        if edgeSize2*2 > sqrt(2)*HMax
            failure = 0;
        else
            edgeSize1 = edgeSize2;
            edgeSize2 = edgeSize1*2;
        end
    else
        failure = 0;
    end
end
edgeSize = edgeSize1;

% Build p
Xcoord = zeros(SizeX/edgeSize,2*SizeY/edgeSize);
Ycoord = zeros(SizeX/edgeSize,2*SizeY/edgeSize);
for i = 1:SizeX/edgeSize+1
    for j = 1:2*SizeY/edgeSize+1
        Xcoord(i,j) = 0+(i-1)*edgeSize;
        Ycoord(i,j) = -SizeY+(j-1)*edgeSize;
    end
end
p = [Xcoord(:),Ycoord(:)]';

% Build t
t = zeros(3,2*SizeX/edgeSize*2*SizeY/edgeSize);
k = 1;
for j = 1:2*SizeY/edgeSize
    for i = 1:SizeX/edgeSize
        t(1,k) = i+(j-1)*(SizeX/edgeSize+1);
        t(2,k) = i+1+(j-1)*(SizeX/edgeSize+1);
        t(3,k) = i+j*(SizeX/edgeSize+1);
        k = k+1;
        
        t(1,k) = i+1+(j-1)*(SizeX/edgeSize+1);
        t(2,k) = i+1+j*(SizeX/edgeSize+1);
        t(3,k) = i+j*(SizeX/edgeSize+1);
        k = k+1;
    end
end


%% Identify fully-resolved regions
% UNIFORM LOADING TEST
if strcmp(ExampleType,'uniform_loading')
    [Xa,Ya] = meshgrid(0:dSizeX:FineX,-FineY:dSizeY:FineY);
    
    % PURE BENDING TEST
elseif strcmp(ExampleType,'pure_bending')
    [Xa,Ya] = meshgrid(0:dSizeX:FineX,-SizeY:dSizeY:-SizeY+FineY);
    
    % INDENTATION
elseif strcmp(ExampleType,'indentation')
    [Xa,Ya] = meshgrid(0:dSizeX:FineX,SizeY-FineY:dSizeY:SizeY);
end

%% Refine the initial regular mesh
% Until atoms_refine are repatoms
fail = 1;
while fail
    triangles_to_refine = [];
    for i = 1:size(t,2)
        P1 = p(:,t(1,i));
        P2 = p(:,t(2,i));
        P3 = p(:,t(3,i));
        if norm(cross([P3-P1;0],[P2-P1;0])) > dSizeX^2+TOL_g
            % Test if any atom from IDPad is within current triangle
            Bar1 = ((P2(2) - P3(2))*(Xa - P3(1)) + (P3(1) - P2(1))*...
                (Ya - P3(2))) ./ ((P2(2) - P3(2))*(P1(1) - P3(1)) +...
                (P3(1) - P2(1))*(P1(2) - P3(2)));
            Bar2 = ((P3(2) - P1(2))*(Xa - P3(1)) + (P1(1) - P3(1))*...
                (Ya - P3(2))) ./ ((P2(2) - P3(2))*(P1(1) - P3(1)) +...
                (P3(1) - P2(1))*(P1(2) - P3(2)));
            Bar3 = 1 - Bar1 - Bar2;
            % If yes, add it to triangles_to_refine
            if find(Bar1>-TOL_g & Bar2>TOL_g & Bar3>TOL_g |...
                    Bar1>TOL_g & Bar2>-TOL_g & Bar3>TOL_g |...
                    Bar1>TOL_g & Bar2>TOL_g & Bar3>-TOL_g,1)
                triangles_to_refine = union(triangles_to_refine,i);
            end
        end
    end
    
    % Refine
    [t,p,~] = regular_refine(p,t(1:3,:),triangles_to_refine,dSizeX,TOL_g);
    
    % Try to find all Xa, Ya atoms in p
    fail = 0;
    for i = 1:length(Xa)
        if min(abs(p(1,:)-Xa(i))+abs(p(2,:)-Ya(i)))>TOL_g
            fail = 1;
            break;
        end
    end
end

%% Take mirror image with respect to y-axis
Mt = t(1:3,:)+size(p,2);
Mp = [-p(1,:)
    p(2,:)];
t = [t(1:3,:),Mt];
p = [p,Mp];
fprintf('time consumed %g s\n',toc);

end
