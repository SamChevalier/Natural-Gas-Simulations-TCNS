function [SPFCR] = GAS_Build_SPFCR(Bus,Line)
% GAS_Build_SPFCR:    Build the slack presure function, as a function of
%                     compressor ratio values
%
% Input:
% 1) Bus       Bus structure
% 2) Line      Line structure
%
% Output:
% 1) SPFCR     Numerical value of the source pressure (function)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bus and line counts
m = size(Bus.Inc_M,1); % # of lines
n = size(Bus.Inc_M,2); % # of buses

% Build symbolic vector
CR = sym('CR',[n 1]);

% Build necessary matrix structures
K0 =  0.5*(abs(Bus.Inc_M)' + (Bus.Inc_M)');
KL = -0.5*(abs(Bus.Inc_M)' - (Bus.Inc_M)');

% Define Gamma matrices
Gam1 = diag(Line.Length/2)*(K0'*diag(CR) - KL');
Gam2 = diag(Line.Length/2);
Gam3 = diag(Line.a)^2*(KL'+K0'*diag(CR));
Gam4 = diag(Line.Length.*Line.Lambda./(4*Line.Diam));
Gam5 = K0'*diag(CR) - KL';

% Add to structure
Gam.G3 = Gam3;
Gam.G4 = Gam4;
Gam.G5 = Gam5;

% Define ODE matrix Mx
Mx  = [zeros(n,n)  K0          KL;
       Gam1        zeros(m,m)  zeros(m,m);
       zeros(m,n)  Gam2        Gam2];

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
p_funcs = matlabFunction(SPs,'vars',{phi_nsd,phi_sd,rho_ns,phi_s,phi_ns,CR});
phi     = SS(1:(2*m));
rho     = SS((2*m+1):(2*m+n));

% Call particular values
phi_nsdN = 0;
phi_sdN  = 0;
rho_nsN  = rho(nsb);
phi_sN   = phi(Bus.Ln_slk + (m-Bus.sc));  % These two can actually be interchanged
phi_nsN  = phi(Bus.Ln_slk + (0+Bus.sc));
rho_SS   = SS((2*m+1):(2*m+n));
p_vals   = p_funcs(phi_nsdN,phi_sdN,rho_nsN,phi_sN,phi_nsN,Bus.Comp_Rat);
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
SPFCR  = matlabFunction(SP,'vars',{phi_nsd,phi_sd,rhoV,phi_0V,phi_LV,CR});

end

% Local function: f
function vec = f(phi_L,phi_0,Gam5_rho)
    % Mass flux conservation
    vec = (phi_L+phi_0).*abs(phi_L+phi_0)./Gam5_rho;
end