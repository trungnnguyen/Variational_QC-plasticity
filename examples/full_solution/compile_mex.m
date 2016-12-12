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
    
    disp('build_grad_r');
    eval('mex build_grad_r.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('build_hess_r');
    eval('mex build_hess_r.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('return_mapping');
    eval('mex return_mapping.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('build_en');
    eval('mex build_en.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('build_diss');
    eval('mex build_diss.cpp COMPFLAGS="/openmp $COMPFLAGS"');

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
    
    disp('build_grad_r.cpp');
    eval('mex -v -g build_grad_r.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('build_hess_r.cpp');
    eval('mex -v -g build_hess_r.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('return_mapping.cpp');
    eval('mex -v -g return_mapping.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('build_en');
    eval('mex -v -g build_en.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    disp('build_diss');
    eval('mex -v -g build_diss.cpp COMPFLAGS="/openmp $COMPFLAGS"');
    
    % Drawing tool is available only for pc
    if ispc
        disp('call_OpenGl');
        eval('mex -v -g call_OpenGl.cpp');
    end
end
clear debugMode;
cd ..;
