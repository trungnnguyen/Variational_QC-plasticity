%% Compile all the *.mex files
% 0 - compile for release, 1 - compile for debugging
debugMode = 0;
clc;

%% Compiling the source codes
disp('Compiling source code:');
disp('----------------------');
cd mex;
if debugMode == 0
    disp('build_atoms');
    eval('mex build_atoms.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('build_bonds');
    eval('mex build_bonds.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('sort_atoms_QC');
    eval('mex sort_atoms_QC.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('sort_sampling_atoms_QC');
    eval('mex sort_sampling_atoms_QC.cpp');
    
    disp('build_grad_r_QC');
    eval('mex build_grad_r_QC.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('build_hess_r_QC');
    eval('mex build_hess_r_QC.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('return_mapping_QC');
    eval('mex return_mapping_QC.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('build_en_QC');
    eval('mex build_en_QC.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('build_diss_QC');
    eval('mex build_diss_QC.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('refine_mesh');
    eval('mex refine_mesh.cpp');    
    
    % Drawing tool is available only for pc
    if ispc
        disp('call_OpenGl');
        eval('mex call_OpenGl.cpp');
    end
elseif debugMode == 1
    disp('build_atoms.cpp');
    eval('mex -v -g build_atoms.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('build_bonds.cpp');
    eval('mex -v -g build_bonds.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('sort_atoms_QC');
    eval('mex -v -g sort_atoms_QC.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('sort_sampling_atoms_QC');
    eval('mex -v -g sort_sampling_atoms_QC.cpp');
    
    disp('build_grad_r_QC');
    eval('mex -v -g build_grad_r_QC.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('build_hess_r_QC');
    eval('mex -v -g build_hess_r_QC.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('return_mapping_QC');
    eval('mex -v -g return_mapping_QC.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('build_en_QC');
    eval('mex -v -g build_en_QC.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('build_diss_QC');
    eval('mex -v -g build_diss_QC.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('refine_mesh');
    eval('mex -v -g refine_mesh.cpp');s    

    % Drawing tool is available only for pc
    if ispc
        disp('call_OpenGl');
        eval('mex -v -g call_OpenGl.cpp');
    end
end
clear debugMode;
cd ..;
