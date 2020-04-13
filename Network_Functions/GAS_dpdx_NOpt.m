function [dpdx] = GAS_dpdx_NOpt(rho_ns,phi_0,phi_L,Gam,Inputs,Mats,Bus,SPFCR)
% GAS_Num_SPFd:  Generate numerical derivative for slack pressure - this is
%                the slow way to complete this task, because compressor
%                ratios are being continuously altered (NOpt = Non-Optimal)
%
% Input:
% 1) rho_ns    Pressure state variables (not including slack)
% 2) phi_0     Flux (0) state variables
% 3) phi_L     Flux (0) state variables
% 4) SPF       Slack pressure function
% 5) NSFD      Non-slack flux derivative function
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
m  = size(Bus.Inc_M,1);
n  = size(Bus.Inc_M,2);

% Compute base
phi_nsd = GAS_phi_nsd(rho_ns,phi_0,phi_L,Gam,Inputs,Mats.Mx3,Mats.Mx4i,Bus);
phi_sd  = GAS_phi_sd(rho_ns,phi_0,phi_L,Gam,Inputs,Mats.Mx3,Mats.Mx4i,Bus);
rho_0   = SPFCR(phi_nsd,phi_sd,rho_ns,phi_0,phi_L,Gam.CR);

% Define x
x = [rho_ns;phi_0;phi_L];
rp = zeros(1,length(x));
for ii = 1:length(x)
    xp     = x;
    xp(ii) = xp(ii) + ep;
    r      = xp(1:(n-1));
    p0     = xp(n:(n+m-1));
    pL     = xp((n+m):(n+2*m-1));
    
    % Compute perturbed rho
    phi_nsd = GAS_phi_nsd(r,p0,pL,Gam,Inputs,Mats.Mx3,Mats.Mx4i,Bus);
    phi_sd  = GAS_phi_sd(r,p0,pL,Gam,Inputs,Mats.Mx3,Mats.Mx4i,Bus);
    rho_p   = SPFCR(phi_nsd,phi_sd,r,p0,pL,Gam.CR);
    rp(ii)  = rho_p;
end

% Compute perturbation
dpdx = (rp-rho_0)/ep;
end

