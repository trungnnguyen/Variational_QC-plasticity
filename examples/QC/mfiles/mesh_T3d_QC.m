function [p,t] = mesh_T3d_QC(ExampleType,SizeX,SizeY,dSizeX,dSizeY,...
    FineX,FineY)

%% Generate a T3d.in file
fprintf('Generating mesh through T3d...\n');
tic;
fid = fopen('T3d.in','w');

% UNIFORM LOADING TEST
if strcmp(ExampleType,'uniform_loading')
    
    % Triangulate just the upper right rectangle
    % Add coarse region points
    fprintf(fid,'# Rough geometry vertices\n\n');
    fprintf(fid,'vertex 1 xyz %.14e %.14e %.14e\n',FineX,0,0);
    fprintf(fid,'vertex 2 xyz %.14e %.14e %.14e\n',SizeX,0,0);
    fprintf(fid,'vertex 3 xyz %.14e %.14e %.14e\n',SizeX,SizeY,0);
    fprintf(fid,'vertex 4 xyz %.14e %.14e %.14e\n',0,SizeY,0);
    fprintf(fid,'vertex 5 xyz %.14e %.14e %.14e\n',0,FineY,0);
    fprintf(fid,'vertex 6 xyz %.14e %.14e %.14e\n\n',FineX,FineY,0);
    
    % Boundary lines
    fprintf(fid,'# Rough geometry curves\n\n');
    fprintf(fid,'curve 1 vertex 1 2 size %f \n',SizeX);
    fprintf(fid,'curve 2 vertex 2 3 size %f \n',SizeY);
    fprintf(fid,'curve 3 vertex 3 4 size %f \n',SizeY);
    fprintf(fid,'curve 4 vertex 4 5 size %f \n',SizeX);
    fprintf(fid,'curve 5 vertex 5 6 size %f \n',dSizeX);
    fprintf(fid,'curve 6 vertex 6 1 size %f \n',dSizeY);
    
    % Patch
    fprintf(fid,'# Rough geometry patches\n\n');
    fprintf(fid,'patch 1 normal 0 0 1 boundary curve 1 2 3 4 5 6\n\n');
    
    % Add all the rectangles in the fully resolved region
    fprintf(fid,'# Fine geometry vertices, lines, surfaces\n\n');
    countV = 10;
    countS = 10;
    countC = 10;
    for i = 1:FineX/dSizeX
        for j = 1:FineY/dSizeY
            
            % Vertices
            fprintf(fid,'vertex %u xyz %.14e %.14e %.14e\n',4*countV-3,...
                (i-1)*dSizeX,(j-1)*dSizeY,0);
            fprintf(fid,'vertex %u xyz %.14e %.14e %.14e\n',4*countV-2,...
                (i)*dSizeX,(j-1)*dSizeY,0);
            fprintf(fid,'vertex %u xyz %.14e %.14e %.14e\n',4*countV-1,...
                (i-1)*dSizeX,(j)*dSizeY,0);
            fprintf(fid,'vertex %u xyz %.14e %.14e %.14e\n\n',4*countV,...
                (i)*dSizeX,(j)*dSizeY,0);
            
            % Curves
            if j == FineY/dSizeY
                fprintf(fid,'curve %u vertex %u %u\n',4*countC-3,...
                    4*countV-3,4*countV-2);
                fprintf(fid,'curve %u vertex %u %u\n',4*countC-2,...
                    4*countV-2,4*countV);
                fprintf(fid,'curve %u vertex %u %u coincide curve 3\n',...
                    4*countC-1,4*countV,4*countV-1);
                fprintf(fid,'curve %u vertex %u %u\n\n',4*countC,...
                    4*countV-1,4*countV-3);
            elseif i == FineX/dSizeX
                fprintf(fid,'curve %u vertex %u %u\n',4*countC-3,...
                    4*countV-3,4*countV-2);
                fprintf(fid,'curve %u vertex %u %u coincide curve 2\n',...
                    4*countC-2,4*countV-2,4*countV);
                fprintf(fid,'curve %u vertex %u %u\n',4*countC-1,...
                    4*countV,4*countV-1);
                fprintf(fid,'curve %u vertex %u %u\n\n',4*countC,...
                    4*countV-1,4*countV-3);
            else
                fprintf(fid,'curve %u vertex %u %u\n',4*countC-3,...
                    4*countV-3,4*countV-2);
                fprintf(fid,'curve %u vertex %u %u\n',4*countC-2,...
                    4*countV-2,4*countV);
                fprintf(fid,'curve %u vertex %u %u\n',4*countC-1,...
                    4*countV,4*countV-1);
                fprintf(fid,'curve %u vertex %u %u\n\n',4*countC,...
                    4*countV-1,4*countV-3);
            end
            
            % Surfaces
            fprintf(fid,'surface %u curve %u %u %u %u\n\n',countS,...
                4*countC-3,4*countC-2,4*countC-1,4*countC);
            
            % Count
            countV = countV + 1;
            countS = countS + 1;
            countC = countC + 1;
        end
    end
    
    % PURE BENDING TEST OR INDENTATION
elseif strcmp(ExampleType,'pure_bending') ||...
        strcmp(ExampleType,'indentation')
    
    % Triangulate just the right rectangle
    % Add coarse region points
    fprintf(fid,'# Rough geometry vertices\n\n');
    fprintf(fid,'vertex 1 xyz %.14e %.14e %.14e\n',FineX,-SizeY,0);
    fprintf(fid,'vertex 2 xyz %.14e %.14e %.14e\n',SizeX,-SizeY,0);
    fprintf(fid,'vertex 3 xyz %.14e %.14e %.14e\n',SizeX,SizeY,0);
    fprintf(fid,'vertex 4 xyz %.14e %.14e %.14e\n',0,SizeY,0);
    fprintf(fid,'vertex 5 xyz %.14e %.14e %.14e\n',0,-SizeY+FineY,0);
    fprintf(fid,'vertex 6 xyz %.14e %.14e %.14e\n\n',FineX,-SizeY+FineY,0);
    
    % Boundary lines
    fprintf(fid,'# Rough geometry curves\n\n');
    fprintf(fid,'curve 1 vertex 1 2 size %f \n',SizeX);
    fprintf(fid,'curve 2 vertex 2 3 size %f \n',SizeY);
    fprintf(fid,'curve 3 vertex 3 4 size %f \n',SizeY);
    fprintf(fid,'curve 4 vertex 4 5 size %f \n',SizeX);
    fprintf(fid,'curve 5 vertex 5 6 size %f \n',dSizeX);
    fprintf(fid,'curve 6 vertex 6 1 size %f \n',dSizeY);
    
    % Patch
    fprintf(fid,'# Rough geometry patches\n\n');
    fprintf(fid,'patch 1 normal 0 0 1 boundary curve 1 2 3 4 5 6\n\n');
    
    % Add all the rectangles in the fully resolved region
    fprintf(fid,'# Fine geometry vertices, lines, surfaces\n\n');
    countV = 10;
    countS = 10;
    countC = 10;
    for i = 1:FineX/dSizeX
        for j = 1:FineY/dSizeY
            
            % Vertices
            fprintf(fid,'vertex %u xyz %.14e %.14e %.14e\n',4*countV-3,...
                (i-1)*dSizeX,(j-1)*dSizeY-SizeY,0);
            fprintf(fid,'vertex %u xyz %.14e %.14e %.14e\n',4*countV-2,...
                (i)*dSizeX,(j-1)*dSizeY-SizeY,0);
            fprintf(fid,'vertex %u xyz %.14e %.14e %.14e\n',4*countV-1,...
                (i-1)*dSizeX,(j)*dSizeY-SizeY,0);
            fprintf(fid,'vertex %u xyz %.14e %.14e %.14e\n\n',4*countV,...
                (i)*dSizeX,(j)*dSizeY-SizeY,0);
            
            % Curves
            if j == FineY/dSizeY
                fprintf(fid,'curve %u vertex %u %u\n',4*countC-3,...
                    4*countV-3,4*countV-2);
                fprintf(fid,'curve %u vertex %u %u\n',4*countC-2,...
                    4*countV-2,4*countV);
                fprintf(fid,'curve %u vertex %u %u coincide curve 3\n',...
                    4*countC-1,4*countV,4*countV-1);
                fprintf(fid,'curve %u vertex %u %u\n\n',4*countC,...
                    4*countV-1,4*countV-3);
            elseif i == FineX/dSizeX
                fprintf(fid,'curve %u vertex %u %u\n',4*countC-3,...
                    4*countV-3,4*countV-2);
                fprintf(fid,'curve %u vertex %u %u coincide curve 2\n',...
                    4*countC-2,4*countV-2,4*countV);
                fprintf(fid,'curve %u vertex %u %u\n',4*countC-1,...
                    4*countV,4*countV-1);
                fprintf(fid,'curve %u vertex %u %u\n\n',4*countC,...
                    4*countV-1,4*countV-3);
            else
                fprintf(fid,'curve %u vertex %u %u\n',4*countC-3,...
                    4*countV-3,4*countV-2);
                fprintf(fid,'curve %u vertex %u %u\n',4*countC-2,...
                    4*countV-2,4*countV);
                fprintf(fid,'curve %u vertex %u %u\n',4*countC-1,...
                    4*countV,4*countV-1);
                fprintf(fid,'curve %u vertex %u %u\n\n',4*countC,...
                    4*countV-1,4*countV-3);
            end
            
            % Surfaces
            fprintf(fid,'surface %u curve %u %u %u %u\n\n',countS,...
                4*countC-3,4*countC-2,4*countC-1,4*countC);
            
            % Count
            countV = countV + 1;
            countS = countS + 1;
            countC = countC + 1;
        end
    end
end
fclose('all');

%% T3d call
dos(['T3d.exe -i T3d.in -o T3d.out -E 1.5 -d ',num2str(max(SizeX,SizeY))]);

%% Read T3d.out
fid = fopen('T3d.out','r');
a = fscanf(fid,'%u',[1 4]);
b = fscanf(fid,'%u',[1 8]);
c = fscanf(fid,'%f',[a(1), b(1)])';
d = fscanf(fid,'%f',[6, b(2)])';
e = fscanf(fid,'%f',[a(1), b(3)])';
fclose('all');

% Construct p and t matrices
p = c(:,2:3)';
t = e(:,2:4)';

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
