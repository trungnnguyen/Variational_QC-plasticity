function [t,p,meshConverged] = regular_refine(p,t,refine_triangles,...
    dSize,TOL_g)

% Check for fully resolved triangles
indDelete = [];
for i = 1:length(refine_triangles)
    P1 = [p(1,t(1,refine_triangles(i))),p(2,t(1,refine_triangles(i)))]';
    P2 = [p(1,t(2,refine_triangles(i))),p(2,t(2,refine_triangles(i)))]';
    P3 = [p(1,t(3,refine_triangles(i))),p(2,t(3,refine_triangles(i)))]';
    
    darea = norm(cross([P3-P1;0],[P2-P1;0]));
    if darea <= dSize^2+TOL_g
        indDelete = union(indDelete,i);
    end
end
refine_triangles(indDelete) = [];

% Refine the mesh
if ~isempty(refine_triangles)
    [p,t] = refine_mesh(p,t(1:3,:),refine_triangles,1,TOL_g);
    meshConverged = 0;
else
    meshConverged = 1;
end

end
