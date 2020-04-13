function [p_s] = GAS_SPn(Gam,Inputs,Mats,Bus,rho_ns,phi_0,phi_L)
% GAS_SPn:    Build the source pressure function - solve numerically.
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
%                   etc..
% 7) SS        Steady state initial condition structure
%                   phi0       => Initial mass fluxes (left side of line)
%                   phiL       => Initial mass fluxes (right side of line)
%                   rho        => Initial pressures
%                   d_inj      => Initial load/source injections
%
% Output:
% 1) p_s       Numerical value of the source pressure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bus and line counts
m = size(Bus.Inc_M,1); % # of lines
n = size(Bus.Inc_M,2); % # of buses

% Define a symbolic vector, and fill with n-1 numerical value and 1 symbol
syms rho_s
rho             = sym('a',  [n 1]);
v               = 1:n;
v(Bus.Slack)    = [];
rho(v)          = rho_ns;
rho(Bus.Slack)  = rho_s;

% Define full vector of functions
Gt1 = Inputs.MFd_CF;
Gt2 = phi_0-phi_L;
Gt3 = Gam.G3*rho - Gam.G4*f(phi_L,phi_0,Gam.G5*rho);

% Remove the slack bus line momentum equation
Gt3(Bus.Ln_slk) = [];
Gt              = [Gt1; Gt2; Gt3];

% Get the expression for source flux derivative: (Mx4i*(Gt - Mx3*rho_sd)), 
% but (Mx4i*Gt) should be sufficient because phi_0d doesn't depend on the 
% source pressure derivative, so (Mx3*p_sd) = 0 at phi_0d
M      = Mats.Mx4i*Gt;
phi_0d = M(n-1 + m + Bus.Ln_slk);  % 'n-1' because source pressure removed

% Get the element associated with the outgoing flux derivative at the
% source bus (left side flux)
Mx2e = Mats.Mx(n+m+Bus.Ln_slk,n+Bus.Ln_slk);

% Now, build an expression with the source pressure
fv  = (Gam.G3*rho - Gam.G4*f(phi_L,phi_0,Gam.G5*rho));
eq  =   0 == Mx2e*phi_0d - fv(Bus.Ln_slk);
rho_0   = Bus.SS((2*m+1):(2*m+n));
p_s     = vpasolve(eq,rho_s,rho_0(Bus.Slack));

% Choose the solution closest to the steady state value
rho_0   = Bus.SS((2*m+1):(2*m+n));
p_diffs = p_s - rho_0(Bus.Slack);
[~,ind] = min(abs(p_diffs));
p_s     = double(p_s(ind));

end

% Local function: f
function vec = f(phi_L,phi_0,Gam5_rho)
    % Mass flux conservation
    vec = (phi_L+phi_0).*abs(phi_L+phi_0)./Gam5_rho;
end