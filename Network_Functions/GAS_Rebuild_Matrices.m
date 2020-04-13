function [Gam,Mats] = GAS_Rebuild_Matrices(Bus,Line,t)
% GAS_Rebuild_Matrices:  Rebuild system matrices as compressor values are
%                        updated, as parameterized by time "t".
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

% Update
CR     = Bus.Comp_Rat + Bus.CR_factors.*CR_Update(t,Bus.CR_ts,Bus.CR_a);
Gam.CR = CR;

% To test this (for compressor ratio at bus 3 updated), try the following:
% % t = 0:0.1:100;
% % plot(t,Bus.Comp_Rat(3) + Bus.CR_factors(3).*CR_Update(t,Bus.CR_ts,Bus.CR_a));

% Bus and line counts
m = size(Bus.Inc_M,1); % # of lines
n = size(Bus.Inc_M,2); % # of buses

% Find the line tied to the slack bus
Bus.Ln_slk = find(Bus.Inc_M(:,Bus.Slack));

% Test the value of the connection (sc = slack connection)
if Bus.Inc_M(Bus.Ln_slk,Bus.Slack) == 1
    Bus.sc = m;
elseif Bus.Inc_M(Bus.Ln_slk,Bus.Slack) == -1
    Bus.sc = 0;
end

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
Mats.Mx  = [zeros(n,n)  K0          KL;
            Gam1        zeros(m,m)  zeros(m,m);
            zeros(m,n)  Gam2        Gam2];

if strcmp(Bus.Slack_Type,'CP') || isempty(Bus.Slack_Type)  % Option 1: Constant Pressure
    % Segment Mx: (1) remove the equations associated with the conservation law
    % at the specified slack bus
    Mats.Mx4              = Mats.Mx;
    Mats.Mx4(Bus.Slack,:) = [];

    % Segment Mx: (2) remove the dependence on source pressure (column), but
    % save these to subtract to the RHS
    Mats.Mx3              = Mats.Mx4(:,Bus.Slack);
    Mats.Mx4(:,Bus.Slack) = [];
    Mats.Mx4i             = inv(Mats.Mx4);
    
elseif strcmp(Bus.Slack_Type,'SF') % Option 2: Sigmoid Function
    % Segment Mx: (1) remove the equations associated with the conservation law
    % at the specified slack bus
    Mats.Mx4              = Mats.Mx;
    Mats.Mx4(Bus.Slack,:) = [];

    % Segment Mx: (2) remove the dependence on source pressure (column)
    Mats.Mx3              = Mats.Mx4(:,Bus.Slack);
    Mats.Mx4(:,Bus.Slack) = [];
    
elseif strcmp(Bus.Slack_Type,'CF') % Option 3: Constant Flux

    % Segment Mx: (1) remove the equations associated with the momentum
    % at the specified slack bus (n conservation, m mass, and m momentum)
    Mats.Mx4                   = Mats.Mx;
    Mats.Mx4(n+m+Bus.Ln_slk,:) = [];

    % Segment Mx: (2) remove the dependence on source pressure (column), but
    % save these to subtract to the RHS (n pressure, m flux left, m flux
    % right)
    Mats.Mx3              = Mats.Mx4(:,Bus.Slack);
    Mats.Mx4(:,Bus.Slack) = [];
    Mats.Mx4i             = inv(Mats.Mx4);
    
end
end

% Local function: f
function CR = CR_Update(t,CR_ts,CR_a)
    % Mass flux conservation
    CR = exp((t-CR_ts)*CR_a)./(1 + exp((t-CR_ts)*CR_a));
end