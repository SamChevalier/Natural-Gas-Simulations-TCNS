function [SFD,NSFD] = GAS_Build_FD(Gam,Inputs,Mx3,Mx4i,Bus)
% GAS_Build_SPFn:   Build functions to compute slack and non-slack bus flux
%                   derivative values.
%
% Input:
% 1) Gam       Gamma matrices structure
%                   G3  => Gamma3
%                   G4  => Gamma4
%                   G5  => Gamma5
% 2) Inputs    ODE model inputs structure
%                   MFd_CF  => Vector of (n) mass flux injection derivatives
%                   etc..
% 4) Mx3       Portion of ODE coefficient matrix associated with the source
%              source pressure derivative
% 5) Mx4i      Inverted portion of ODE coefficient matrix
% 6) Bus       Node structure with the following elements
%                   Inc_M      => Directed incidence matrix (MxN)
%                   Comp_Rat   => Compressor ratios
%                   Types      => 1=source, 2=load
%                   Slack      => Index of source node which should be the
%                                 constant pressure node, and where we will
%                                 eliminate the conservation law
%                   Slack_Type => The slack source may be constant pressure 
%                                 ('CP') and variable flux (default), or 
%                                 constant flux ('CF') and variable 
%                                 pressure.
%                   etc..
%
% Output:
% 1) SFD       Slack flux derivative
% 2) NSFD      Non-slack flux derivative
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bus and line counts
m = size(Bus.Inc_M,1); % # of lines
n = size(Bus.Inc_M,2); % # of buses

% Create vectors of symbolic vectors
rho_ns = sym('rho_ns', [n-1 1]);
phi_0  = sym('phi_0',  [m 1]);
phi_L  = sym('phi_L',  [m 1]);

% Leverage the function GAS_phi_nsd
[phi_sd]  = GAS_phi_sd(rho_ns,phi_0,phi_L,Gam,Inputs,Mx3,Mx4i,Bus);
[phi_nsd] = GAS_phi_nsd(rho_ns,phi_0,phi_L,Gam,Inputs,Mx3,Mx4i,Bus);

% Create output function
SFD   = matlabFunction(phi_sd,'vars',{rho_ns,phi_0,phi_L});
NSFD  = matlabFunction(phi_nsd,'vars',{rho_ns,phi_0,phi_L});

end
