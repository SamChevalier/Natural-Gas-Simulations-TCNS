function [r,Bus] = GAS_Dyn_Sim(IC,Bus,Line,Inputs,TD_Sim)
% GAS_DYN_SIM:    Simulate a Gas Network. In this framework, we assume there
%                 are (n+2m) state variables: two mass flux variables on each
%                 line and 1 pressure variable associated with each node.
%
% Inputs:
% 1) IC        Initial condition structure
%                   phi0       => Initial mass fluxes (left side of line)
%                   phiL       => Initial mass fluxes (right side of line)
%                   rho        => Initial pressures
%                   d_inj      => Initial load/source injections
% 2) Bus       Node structure with the following elements
%                   SS         => Steady State vector which has the same
%                                 structure as the IC structure
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
% 3) Line      Line structure with the following elements
%                   Length     => Directed incidence matrix
%                   Lambda     => 1=source, 2=load
%                   Diam       => Load (2) injections 
%                   a          => Propagation/speed factor
% 4) Inputs    ODE model inputs structure
%                   MFd_CP     => Vector of (n) mass flux injection derivatives
%                   MFd_SF     => Vector of (n) mass flux injection derivatives
%                   MFd_CF     => Vector of (n) mass flux injection derivatives
% 5) TD_Sim    Simulation characteristics structure
%                   dt         => Time step
%                   tf         => Final time
%
% Outputs:
% 1) r         Simulation results structure
%                   t          => Output time
%                   phi0       => Mass flow at x = 0
%                   phiL       => Mass flow at x = L
%                   rho        => Nodal pressure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bus and line counts
m = size(Bus.Inc_M,1); % # of lines
n = size(Bus.Inc_M,2); % # of buses

% Find the line tied to the slack bus
Bus.Ln_slk = find(Bus.Inc_M(:,Bus.Slack));
% If the slack node is tied to two (or more) branches, arbitrarily use the
% first branch. All branches will represent valid ways to compute the
% slack pressure. If the slack isn't a terminal node, then its flux
% derivative will be nonzero even though the flux injection is constant.
if length(Bus.Ln_slk)>1
    Bus.Ln_slk = Bus.Ln_slk(1);
end

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
Gam1 = diag(Line.Length/2)*(K0'*diag(Bus.Comp_Rat) - KL');
Gam2 = diag(Line.Length/2);
Gam3 = diag(Line.a)^2*(KL'+K0'*diag(Bus.Comp_Rat));
Gam4 = diag(Line.Length.*Line.Lambda./(4*Line.Diam));
Gam5 = K0'*diag(Bus.Comp_Rat) - KL';

% Add to structure
Gam.G3 = Gam3;
Gam.G4 = Gam4;
Gam.G5 = Gam5;

% Define ODE matrix Mx
Mats.Mx  = [zeros(n,n)  K0          KL;
            Gam1        zeros(m,m)  zeros(m,m);
            zeros(m,n)  Gam2        Gam2];

% If specified, hold flux at the source bus constant and let pressure vary
if strcmp(Bus.Slack_Type,'CP') || isempty(Bus.Slack_Type)  % Option 1: Constant Pressure

    % Set up Indices
    nslk_inds            = 1:n;
    nslk_inds(Bus.Slack) = [];

    % Segment Mx: (1) remove the equations associated with the conservation law
    % at the specified slack bus
    Mats.Mx4              = Mats.Mx;
    Mats.Mx4(Bus.Slack,:) = [];

    % Segment Mx: (2) remove the dependence on source pressure (column), but
    % save these to subtract to the RHS
    Mats.Mx3              = Mats.Mx4(:,Bus.Slack);
    Mats.Mx4(:,Bus.Slack) = [];
    Mats.Mx4i             = inv(Mats.Mx4);

    % Define initial conditions & time vector
    phi0_0 = IC(1:m);
    phiL_0 = IC((m+1):(2*m));
    rho_0  = IC((2*m+1):(2*m+n));
    y0     = [rho_0(nslk_inds); phi0_0; phiL_0];
    tvec   = 0:TD_Sim.dt:TD_Sim.tf;
    
    % Set slack pressure to its stead-state value (this may be altered if
    % a custom input is desired)
    Inputs.ps  = rho_0(Bus.Slack);   % pressure
    Inputs.psd = 0;                  % pressure derivative
    
    % Populate (unused) Funcs structure
    Funcs = struct('f',0);
    % % opt   = odeset('RelTol',1e-7,'AbsTol',1e-7);

    % Simulate with ode15s (or ode45)
    tic
    [t_out,x_out] = ode23tb(@(t,x) GAS_ODEs(t,x,Inputs,m,n,Mats,Gam,Bus,nslk_inds,Funcs,Line),tvec,y0);
    toc

    % Parse output
    r.t    = t_out;
    r.rho  = x_out(:,1:(n-1));
    r.phi0 = x_out(:,n:(m+n-1));
    r.phiL = x_out(:,(n+m):(2*m+n-1));
    
    % Add the source pressure (although constant) to the vector r.rho
    r.rho = [r.rho(:,1:(Bus.Slack-1)) Inputs.ps*ones(length(t_out),1) r.rho(:,(Bus.Slack):end)];
    
elseif strcmp(Bus.Slack_Type,'SF') % Option 2: Sigmoid Function
    
    % Call parameters
    gamma   = Bus.gamma;
    phi_max = Bus.phi_max;
    phi_MAX = Bus.phi_MAX;
    
    % Compute "nominal" rho value: as currently written, it only make sense
    % to have the sigmoid function applied at a source bus which is also a
    % terminal bus (i.e. starting bus)
    if isfield(Bus, 'rho_nom')
        rho_nom = Bus.rho_nom;
    else
        rhoS_0  = IC(2*m+Bus.Slack);
        phiS_0  = IC(Bus.Ln_slk);
        rho_nom = rhoS_0/(exp(gamma*(phi_MAX - phiS_0))/(1+exp(gamma*(phi_MAX - phiS_0))));
    end
    
    % What is the initial value of z?
    phiS_0 = IC(Bus.Ln_slk);
    z0     = log( (phiS_0/phi_max)/(1-(phiS_0/phi_max)) );
    
    % Build Sigmoid Functions
    [h1,h2,Phi,Rho] = GAS_Build_Sig(gamma,phi_max,rho_nom,phi_MAX);
    
    % Populate Funcs structure
    Funcs.h1  = h1;
    Funcs.h2  = h2;
    Funcs.Phi = Phi;
    Funcs.Rho = Rho;

    % Set up Indices
    nslk_inds            = 1:n;
    nslk_inds(Bus.Slack) = [];

    % Segment Mx: (1) remove the equations associated with the conservation law
    % at the specified slack bus
    Mats.Mx4              = Mats.Mx;
    Mats.Mx4(Bus.Slack,:) = [];

    % Segment Mx: (2) remove the dependence on source pressure (column)
    Mats.Mx3              = Mats.Mx4(:,Bus.Slack);
    Mats.Mx4(:,Bus.Slack) = [];

    % Define initial conditions & time vector
    phi0_0 = IC(1:m);
    phiL_0 = IC((m+1):(2*m));
    rho_0  = IC((2*m+1):(2*m+n));
    y0     = [rho_0(nslk_inds); phi0_0; phiL_0; z0]; % Append zero
    tvec   = 0:TD_Sim.dt:TD_Sim.tf;
    
    % Simulate with ode15s (or ode45 [ode23t can be unstable])
    tic
    [t_out,x_out] = ode23tb(@(t,x) GAS_ODEs(t,x,Inputs,m,n,Mats,Gam,Bus,nslk_inds,Funcs,Line),tvec,y0);
    toc
    
    % Parse output (doesn't include source pressure yet)
    r.t       = t_out;
    r.phi0    = x_out(:,n:(m+n-1));
    r.phiL    = x_out(:,(n+m):(2*m+n-1));
    r.z       = x_out(:,end);
    r.rho_nom = rho_nom;
    
    % Parse pressure (add source pressure)
    phi_slk            = x_out(:,n-1 + Bus.Ln_slk + (m - Bus.sc));
    Rho_S              = Rho(phi_slk);
    r.rho              = zeros(length(r.t),n);
    r.rho(:,nslk_inds) = x_out(:,1:(n-1));
    r.rho(:,Bus.Slack) = Rho_S;

elseif strcmp(Bus.Slack_Type,'CF') % Option 3: Constant Flux
    
    % Set up Indices
    nslk_inds            = 1:n;
    nslk_inds(Bus.Slack) = [];

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

    % Compute the slack bus pressure via analytical solution of QF
    [SPF]       = GAS_Build_SPF(Gam,Mats.Mx,Bus);
    [SFD,NSFD]  = GAS_Build_FD(Gam,Inputs,Mats.Mx3,Mats.Mx4i,Bus);
    
    % A quick note about SPF: it generates two solutions (via solutions to
    % the quadratic formula), and the correct one is found by comparing the
    % function to the "SS" values. If a parameter, such as a compressor
    % ratio, is altered, it is possible that is picks the "wrong" solution.
    % Watch out for this.
    
    % Previous solution to this problem (untractable for large systems)
    % % [SPF]     = GAS_Build_SPFQA(Gam,Mats.Mx,Bus,SS);
    % % [SPdF]    = GAS_Build_SPdF(Gam,Inputs,Mats.Mx,Mats.Mx3,Mats.Mx4i,Bus);
    % % [phi_0dF] = GAS_Build_phi_0dF(Gam,Inputs,Mats.Mx3,Mats.Mx4i,Bus);

    % Populate Funcs structure
    Funcs.SPF     = SPF;
    Funcs.SFD     = SFD;
    Funcs.NSFD    = NSFD;
    
    % Define initial conditions & time vector
    phi0_0 = IC(1:m);
    phiL_0 = IC((m+1):(2*m));
    rho_0  = IC((2*m+1):(2*m+n));
    y0     = [rho_0(nslk_inds); phi0_0; phiL_0];  % Slack pressure removed
    tvec   = 0:TD_Sim.dt:TD_Sim.tf;
    
    % Are compressor ratios being updated?
    if Bus.CR_update == 1
        [SPFCR]     = GAS_Build_SPFCR(Bus,Line);
        Funcs.SPFCR = SPFCR;
    end
    
    % Define a pre-conditioner, if helpful
    % %     dpdx    = GAS_dpdx(rho_0(nslk_inds),phi0_0,phiL_0,SPF,NSFD);
    % %     MM      = [Mats.Mx3  Mats.Mx4; 1  -dpdx];
    % %     Mats.PC = MM\eye(size(MM));

    % Simulate with ode23tb (or ode15s)
    % % % opt = odeset('RelTol',1e-6,'AbsTol',1e-6); % % %
    tic
    [t_out,x_out] = ode23tb(@(t,x) GAS_ODEs(t,x,Inputs,m,n,Mats,Gam,Bus,nslk_inds,Funcs,Line),tvec,y0);
    toc
    
    % Parse output (doesn't include source pressure yet)
    r.t    = t_out;
    r.rho  = x_out(:,1:(n-1));
    r.phi0 = x_out(:,n:(m+n-1));
    r.phiL = x_out(:,(n+m):(2*m+n-1));
    
    % Post-Computation: are compressor ratios being updated?
    if Bus.CR_update == 1
        sp = zeros(length(t_out),1);
        % This is a faily slow routine
        for ii  = 1:length(t_out)
            t          = r.t(ii);
            [Gam,Mats] = GAS_Rebuild_Matrices(Bus,Line,t);
            rns        = r.rho(ii,:)';
            p0         = r.phi0(ii,:)';
            pL         = r.phiL(ii,:)';
            pnsd       = GAS_phi_nsd(rns,p0,pL,Gam,Inputs,Mats.Mx3,Mats.Mx4i,Bus);
            psd        = GAS_phi_sd(rns,p0,pL,Gam,Inputs,Mats.Mx3,Mats.Mx4i,Bus);
            sp(ii)     = SPFCR(pnsd,psd,rns,p0,pL,Gam.CR);
        end
    else
        % Post Compute Slack Pressure
        phi_sd  = SFD(r.rho',r.phi0',r.phiL');
        phi_nsd = NSFD(r.rho',r.phi0',r.phiL');
        if length(phi_sd) == 1
            phi_sd = phi_sd*ones(size(phi_nsd));
        end
        sp = SPF(phi_nsd,phi_sd,r.rho',r.phi0',r.phiL');
    end
    
    % Update pressure structure
    r.rho              = zeros(length(r.t),n);
    r.rho(:,nslk_inds) = x_out(:,1:(n-1));
    r.rho(:,Bus.Slack) = sp;

    % Append steady state values
    r.t     = [-TD_Sim.dt; r.t];
    r.rho   = [rho_0'; r.rho];
    r.phi0  = [phi0_0'; r.phi0];
    r.phiL  = [phiL_0'; r.phiL];
end
end