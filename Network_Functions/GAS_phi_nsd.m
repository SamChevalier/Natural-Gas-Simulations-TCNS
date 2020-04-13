function [phi_nsd] = GAS_phi_nsd(rho_ns,phi_0,phi_L,Gam,Inputs,Mx3,Mx4i,Bus)
% GAS_phi_nsd: Compute non-slack flux derivative
%
% Input:
% 1) Gam       Gamma matrices structure
%                   G3  => Gamma3
%                   G4  => Gamma4
%                   G5  => Gamma5
% 2) Inputs    ODE model inputs structure
%                   MFd_CF  => Vector of (n) mass flux injection derivatives
%                   etc..
% 3) Mx3       Portion of ODE coefficient matrix associated with the source
%              source pressure derivative
% 4) Mx4i      Inverted portion of ODE coefficient matrix
% 5) Bus       Node structure with the following elements
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
% 1) phi_nsd   Function to compute slack line pressure derivative.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bus and line counts
m = size(Bus.Inc_M,1); % # of lines
n = size(Bus.Inc_M,2); % # of buses

% Numerical matrix: n-1 densities, Bus.Ln_slk^th eq, and Bus.sc=m to flip
Mx4iV = Mx4i(n-1 + Bus.sc + Bus.Ln_slk,:);

% Define full vector of functions
Gt1     = Inputs.MFd_CF;
Gt2     = phi_0-phi_L;
nsb_ind = [1:(Bus.Slack-1) (Bus.Slack+1):n]; % No slack bus
nsl_ind = [1:(Bus.Ln_slk-1) (Bus.Ln_slk+1):m]; % No slack line
G3p     = Gam.G3(nsl_ind,nsb_ind);
G4p     = Gam.G4(nsl_ind,nsl_ind);
G5p     = Gam.G5(nsl_ind,nsb_ind);
Gt3p    = G3p*rho_ns - G4p*f(phi_L(nsl_ind),phi_0(nsl_ind),G5p*rho_ns);

% Assemble
Gt = [Gt1; Gt2; Gt3p];

% Get the expression for source flux derivative: (Mx4i*(Gt - Mx3*rho_sd)), 
% but (Mx4i*Gt) should  be sufficient because phi_0d doesn't depend on the 
% source pressure derivative, so (Mx3*rho_sd) = 0 at the "phi_0d" index
rho_sd  = 0;
phi_nsd = Mx4iV*(Gt - Mx3*rho_sd);

end

% Local function: f
function vec = f(phi_L,phi_0,Gam5_rho)
    % Mass flux conservation
    vec = (phi_L+phi_0).*abs(phi_L+phi_0)./Gam5_rho;
end