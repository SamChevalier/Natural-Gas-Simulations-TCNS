function [phi,rho,d_inj] = GAS_NR_Init(Bus,Line,NR)
% GAS_NR_INIT: Intialize a Gas Network via Newton-Raphson
%
% Inputs:
% 1) Bus       Node structure with the following elements
%                   Inc_M      => Directed incidence matrix (MxN)
%                   Types      => 1=source, 2=load
%                   Src_Press  => Source (1) pressures
%                   Load_Injs  => Load (2) injections 
%                   Comp_Rat   => Compressor ratios
%                   L_param    => *Optional parameter codifying the
%                                 value (D*a^2)/(L*lambda) in a vector of
%                                 scalar values
%                   PG         => *Optional pressure guess vector (length n)
%                   DG         => *Optional injection guess vector (length n)
% 2) Line      Line structure with the following elements
%                   Length     => Line length
%                   Lambda     => Friction factor
%                   Diam       => Line Diameter 
%                   a          => Propagation/speed factor
%                   FG         => *Optional flow guess vector (length m)
%                   Area       => *Optional - if this is present, the flow
%                                  variables are actually kg/s
% 3) NR        Newton Raphson specifications structure
%                   max_its    => Maximum number of allowable iterations
%                   epsilon    => Convergence tolerance
%
% Outputs:
% 1) phi       Complete set (m) of line flow variables
% 2) rho       Complete set (n) of pressures (incoming/left side/pre-compressor)
% 3) d_inj     Complete set (n) of nodal injections
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NR Specifications
its = NR.max_its;
ep  = NR.epsilon;

% Bus and line counts
m = size(Bus.Inc_M,1); % # of lines
n = size(Bus.Inc_M,2); % # of buses

% Indices of sources (alpha) and loads (beta)
source_inds = find(Bus.Types == 1);
    alpha   = length(source_inds);
load_inds   = find(Bus.Types == 2);
    beta    = length(load_inds);

% Build necessary matrix structures
K0 =  0.5*(abs(Bus.Inc_M)' + (Bus.Inc_M)');
KL = -0.5*(abs(Bus.Inc_M)' - (Bus.Inc_M)');

% Use line parameters to define matrix M
if isfield(Line,'L_param')
    M = diag(Line.L_param);
else
    % Build M from diam, length, a, and lambda
    M = diag((Line.Diam.*((Line.a).^2))./(Line.Length.*Line.Lambda));
end

% !!! Actually !!! If line area is specified, then M *must* include the
% line area and flux variables are actually flows in kg/s.
if isfield(Line,'Area')
    M = diag((((Line.Area).^2).*Line.Diam.*((Line.a).^2))./(Line.Length.*Line.Lambda));
end

% Define Gamma0 matrix for Jacobian
Gam0 = M*(K0'*diag(Bus.Comp_Rat).^2 + KL');

% Define and initialize vectors: x  = [phi; rho;  d_inj]
%                                xH = [phi; rhoH; d_injH]
x  = zeros(m+2*n,its);
xH = zeros(m+n,its);

% Define Inds
phi_ind = (1:m)';   rho_ind  = ((m+1):(m+n))';    d_inj_ind  = ((m+n+1):(m+2*n))';
                    rhoH_ind = ((m)+(1:beta))';   d_injH_ind = ((m+beta)+(1:alpha))';
% Add knowns
x(m+source_inds,:) = Bus.Src_Press.*ones(alpha,its);
x(m+n+load_inds,:) = Bus.Load_Injs.*ones(beta,its);

% Initialization: ~~ Flows ~~
if isfield(Line,'FG')
    x(phi_ind,1) = Line.FG; else
    x(phi_ind,1) = sum(abs(Bus.Load_Injs))/n; % Dumb Guess
end

% Initialization: ~~ Pressures ~~
if isfield(Bus,'PG')
    x(rho_ind,1)     = Bus.PG; else
    x(m+load_inds,1) = min(Bus.Src_Press);        % Dumb Guess
end

% Initialization: ~~ Source Injection ~~
if isfield(Bus,'DG')
    x(d_inj_ind,1)       = Bus.DG; else
    x(m+n+source_inds,1) = sum(abs(Bus.Load_Injs))/n; % Dumb Guess
end

% Initial xH
xH(phi_ind,1)    = x(phi_ind,1);
xH(rhoH_ind,1)   = x(m+load_inds,1);
xH(d_injH_ind,1) = x(m+n+source_inds,1);

% Initialize loop
NR_run = 1;
ii     = 2;

% Enter NR loop
while NR_run == 1
    
    % Get current values for Jacobian
    phi = x(phi_ind,ii-1);
    rho = x(rho_ind,ii-1);
    d   = x(d_inj_ind,ii-1);
    
    % Assemble Jacobian
    J11 = (Bus.Inc_M)';
    J12 = zeros(n,n);
    J13 = -eye(n,n);
    J21 = -diag(sign(phi))*2*diag(phi);
    J22 = 2*Gam0*diag(rho);
    J23 = zeros(m,n);
    J   = [J11 J12 J13;
           J21 J22 J23];
    
    % Reduce size of J: Remove Cols of fixed varaibles
    Jh = J(:,[phi_ind; (m+load_inds); (m+n+source_inds)]);
    
    % Compute functions A1 & A2
    A1v = A1(Bus.Inc_M,phi,d);
    A2v = A2(Gam0,rho,phi);
    A   = [A1v; A2v];
    
    % Update xH & x
    inc                   = Jh\A;
    xH(:,ii)              = xH(:,ii-1) - inc;
    x(phi_ind,ii)         = xH(phi_ind,ii);
    x(m+load_inds,ii)     = xH(rhoH_ind,ii);
    x(m+n+source_inds,ii) = xH(d_injH_ind,ii);
    
    % Test
    if (ii >= its) || (max(abs(inc)) < ep)
        NR_run = 0;
    end
        
    % Increment
    ii = ii + 1;
end

% Parse Final Results
phi   = x(phi_ind,ii-1);
rho   = x(rho_ind,ii-1);
d_inj = x(d_inj_ind,ii-1);

end

% Local function: A1
function val = A1(E,phi,d)
    % Mass flux conservation
    val = E'*phi-d;
end

% Local function: A2
function val = A2(Gamma,rho,phi)
    % Mass flux conservation
    val = Gamma*rho.^2-diag(sign(phi))*phi.^2;
end