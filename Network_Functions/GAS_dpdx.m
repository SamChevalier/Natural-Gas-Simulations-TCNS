function [dpdx] = GAS_dpdx(rho_ns,phi_0,phi_L,SPF,SFD,NSFD)
% GAS_Num_SPFd:  Generate numerical derivative for slack pressure
%
% Input:
% 1) rho_ns    Pressure state variables (not including slack)
% 2) phi_0     Flux (0) state variables
% 3) phi_L     Flux (0) state variables
% 4) SPF       Slack pressure function
% 5) SFD       Slack flux derivative function
% 6) NSFD      Non-slack flux derivative function
% 
% Output:
% 1) dpdx      Derivative vector for solving for rho_dot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numerical step
ep = 1e-4;

% Sizes
n = length(rho_ns)+1;
m = length(phi_0);

M  = repmat([rho_ns;phi_0;phi_L],1,n+2*m-1);
pM = M + ep*eye(n+2*m-1);

% Parse
r  = pM(1:(n-1),:);
p0 = pM(n:(n+m-1),:);
pL = pM((n+m):(n+2*m-1),:);

% Compute flux derivative and source pressure for pert/non-pert
phi_sd  = SFD([r rho_ns],[p0 phi_0],[pL phi_L]);
phi_nsd = NSFD([r rho_ns],[p0 phi_0],[pL phi_L]);
rp      = SPF(phi_nsd,phi_sd,[r rho_ns],[p0 phi_0],[pL phi_L]);

% Compute perturbation
dpdx = (rp(1:end-1) - rp(end))/ep;
