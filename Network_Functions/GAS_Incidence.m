function E = GAS_Incidence(ngc)
% GAS_INCIDENCE: Build the incidence matrix associated with the system.
%
% Inputs:
% 1) ngc       Natural gas system structure with the following elements
%                   node       => bus data      
%                   source     => source data
%                   pipeline   => line data
%                   compressor => compressor data
%
% Outputs:
% 1) E         Directed incidence matrix (mxn)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% System size
n = size(ngc.node,1);
m = size(ngc.pipeline,1);

% Build Incidence Matrix
E=zeros(m,n);
for i=1:m
    E(i,ngc.pipeline(i,1))=1;
    E(i,ngc.pipeline(i,2))=-1;
end

end