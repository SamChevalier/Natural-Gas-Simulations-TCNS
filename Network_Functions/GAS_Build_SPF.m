function [SPF] = GAS_Build_SPF(Gam,Mx,Bus)
% GAS_Build_SP:    Build the slack presure function
%
% Input:
% 1) Gam       Gamma matrices structure
%                   G3  => Gamma3
%                   G4  => Gamma4
%                   G5  => Gamma5
% 2) Inputs    ODE model inputs structure
%                   MFd_CF  => Vector of (n) mass flux injection derivatives
%                   etc..
% 3) Mx        ODE coefficient matrix
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
%                   SS         => Steady state initial condition structure
%                       phi0  => Initial mass fluxes (left side of line)
%                       phiL  => Initial mass fluxes (right side of line)
%                       rho   => Initial pressures
%                       d_inj => Initial load/source injections
%                   etc..
%
% Output:
% 1) SPF       Numerical value of the source pressure (function)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bus and line counts
m = size(Bus.Inc_M,1); % # of lines
n = size(Bus.Inc_M,2); % # of buses

% Get SS vector
SS = Bus.SS;

% First, define state variables
syms rho_s
syms rho_ns
syms phi_s
syms phi_ns

% Define symbolic slack/non-slack flux derivatives
syms phi_sd phi_nsd

% Get the element associated with the flux derivative at the
% source bus line
Mx2e_ns = Mx(n+m+Bus.Ln_slk,n+  Bus.sc  +Bus.Ln_slk);
Mx2e_s  = Mx(n+m+Bus.Ln_slk,n+(m-Bus.sc)+Bus.Ln_slk);

% Which bus is the slack? Which bus is the slack tied to?
sb         = Bus.Slack;
s_line     = Bus.Inc_M(Bus.Ln_slk,:);
s_line(sb) = 0;
nsb        = find(s_line);

% Now, build an expression with the source pressure
g1 = Gam.G3(Bus.Ln_slk,sb)*rho_s + Gam.G3(Bus.Ln_slk,nsb)*rho_ns;
g2 = Gam.G5(Bus.Ln_slk,sb)*rho_s + Gam.G5(Bus.Ln_slk,nsb)*rho_ns;
g4 = Gam.G4(Bus.Ln_slk,Bus.Ln_slk);
fv = (g1 - g4*f(phi_ns,phi_s,g2));
eq = 0 == Mx2e_ns*phi_nsd + Mx2e_s*phi_sd - fv;

% This is the quadratic solution of the deleted momentum balance equation
SPs = solve(eq,rho_s);

% Determine the correct solution (heuristically)
p_funcs = matlabFunction(SPs,'vars',{phi_nsd,phi_sd,rho_ns,phi_s,phi_ns});
phi     = SS(1:(2*m));
rho     = SS((2*m+1):(2*m+n));

% Call particular values
phi_nsdN = 0;
phi_sdN  = 0;
rho_nsN  = rho(nsb);
phi_sN   = phi(Bus.Ln_slk + (m-Bus.sc));  % These two can actually be interchanged  
phi_nsN  = phi(Bus.Ln_slk + (0+Bus.sc));
rho_SS   = SS((2*m+1):(2*m+n));
p_vals   = p_funcs(phi_nsdN,phi_sdN,rho_nsN,phi_sN,phi_nsN);
p_diffs  = p_vals - rho_SS(Bus.Slack);

% Select the closest solution
[~,ind]  = min(abs(p_diffs));

% True pressure equation
SP = SPs(ind);

% Add all other state variables
rhoV   = sym('rho',  [n 1]);
phi_0V = sym('phi_0',[m 1]);
phi_LV = sym('phi_L',[m 1]);

% Insert
rhoV(nsb) = rho_ns;
rhoV(sb)  = [];
if Bus.sc == m     % This is also interchangable
    phi_0V(Bus.Ln_slk) = phi_s;
    phi_LV(Bus.Ln_slk) = phi_ns;
else
    phi_0V(Bus.Ln_slk) = phi_ns;
    phi_LV(Bus.Ln_slk) = phi_s;
end

% Create output
SPF  = matlabFunction(SP,'vars',{phi_nsd,phi_sd,rhoV,phi_0V,phi_LV});

end

% Local function: f
function vec = f(phi_L,phi_0,Gam5_rho)
    % Mass flux conservation
    vec = (phi_L+phi_0).*abs(phi_L+phi_0)./Gam5_rho;
end