function dxdt = GAS_ODEs(t,x,Inputs,m,n,Mats,Gam,Bus,nslk_inds,Funcs,Line)
% GAS_ODEs:      Evaluate the gas network ODEs. Use numerical methods
%                method to compute the source pressure.
%
% Inputs:
% 1)  t          Time variable
% 2)  x          State vector x_tilde (rho, phi_0, phi_L)
% 3)  Inputs     ODE model inputs structure
%                   MFd_CP    => Vector of (n-1) mass flux injection derivatives
%                   MFd_CF    => Vector of (n)   mass flux injection derivatives
%                   MFd_SF    => Vector of (n)   mass flux injection derivatives
%                   ps        => Pressure at the slack source
%                   psd       => Pressure derivative at the slack source
% 4)  m          Number of lines
% 5)  n          Number of buses
% 5)  Mats       Matrix structure with matrices
%                   Mx3       => ODE coefficient matrix 3
%                   Mx4i      => ODE coefficient matrix 4 inverse    
% 8)  Gam        Gamma matrices structure
%                   G3        => Gamma3
%                   G4        => Gamma4
%                   G5        => Gamma5
% 9)  Bus        Node structure with the following elements
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
%                   Ln_slk     => The line index associated with the slack
%                                 bus (assume there is only 1)
% 10) nslk_inds  Non-slack bus indices
%                   etc..
% 11) Funcs      Structure of necessary functions
% 12) Line       Line structure
%
% Outputs:
% 1) dxdt        State derivatives for simulator
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% If specified, hold flux at the source bus constant and let pressure vary
if strcmp(Bus.Slack_Type,'CP') || isempty(Bus.Slack_Type)
    
    % Are compressor parameters being updated?
    if Bus.CR_update == 1
        [Gam,Mats] = GAS_Rebuild_Matrices(Bus,Line,t);
    end
    
    % Update
    Mx3  = Mats.Mx3;
    Mx4i = Mats.Mx4i;

    % Define state vectors
    rho             = zeros(n,1);
    rho(Bus.Slack)  = Inputs.ps;
    rho(nslk_inds)  = x(1:(n-1));
    phi_0           = x(n:(m+n-1));
    phi_L           = x((m+n):(2*m+n-1));
    
    % Call ddot
    ddot            = Inputs.MFd_CP;
    ddot(Bus.Slack) = [];
    
    % Split Gt (tilde) into three constituent functions
    Gt1  = ddot;
    Gt2  = phi_0-phi_L;
    Gt3  = Gam.G3*rho - Gam.G4*f(phi_L,phi_0,Gam.G5*rho);
    Gt   = [Gt1; Gt2; Gt3];
    
    % ODE system: set derivative output
    dxdt = Mx4i*(Gt - Mx3*Inputs.psd);
    
elseif strcmp(Bus.Slack_Type,'SF')
    
    % Are compressor parameters being updated?
    if Bus.CR_update == 1
        [Gam,Mats] = GAS_Rebuild_Matrices(Bus,Line,t);
    end
    
    % Inside the state variable vector "x", what is the index of the slack
    % bus flux? There are (n-1) pressure states, the slack line is line
    % number (Bus.Ln_slk), and this index depends on if the slack flux is
    % at the start (0) or end (L) of a line (m-Bus.sc)
    xLn_slk = (n-1) + Bus.Ln_slk + (m-Bus.sc);
    
    % Define z, slack flux, and slack density
    z    = x(end);
    phiS = x(xLn_slk);
    rhoS = Funcs.Rho(phiS);
    
    % Construct M5
    M5            = Mats.Mx4;
    M5(:,xLn_slk) = M5(:,xLn_slk) + Mats.Mx3*Funcs.h2(phiS);
    
    % Construct Ms
    ek          = zeros(1,size(M5,2));
    ek(xLn_slk) = 1;
    Ms          = [ M5  zeros(size(M5,1),1);
                   -ek  Funcs.h1(z)];
    
    % Define state vectors
    rho             = zeros(n,1);
    rho(Bus.Slack)  = rhoS;
    rho(nslk_inds)  = x(1:(n-1));
    phi_0           = x(n:(m+n-1));
    phi_L           = x((m+n):(2*m+n-1));
    z               = x(end);
    
    % Call ddot
    ddot            = Inputs.MFd_SF;
    ddot(Bus.Slack) = [];
    
    % Split Gt (tilde) into three constituent functions
    Gt1 = ddot;
    Gt2 = phi_0-phi_L;
    Gt3 = Gam.G3*rho - Gam.G4*f(phi_L,phi_0,Gam.G5*rho);
    Gt  = [Gt1; Gt2; Gt3; 0]; % Append zero
    
    % ODE system: set derivative output
    dxdt = Ms\Gt;
    
    % Apply scaling function to prevent the blow-up of z
    dxdt(end) = dxdt(end) * 1/(1+exp(z^2 - 100));
    
elseif strcmp(Bus.Slack_Type,'CF')
    
    % Non-source pressures
    rho_ns = x(1:(n-1));
    
    % Define other state vectors
    phi_0 = x(n:(m+n-1));
    phi_L = x((m+n):(2*m+n-1));
    
    % Are compressor parameters being updated?
    if Bus.CR_update == 1
        [Gam,Mats] = GAS_Rebuild_Matrices(Bus,Line,t);
        [dpdx]     = GAS_dpdx_NOpt(rho_ns,phi_0,phi_L,Gam,Inputs,Mats,Bus,Funcs.SPFCR);
    else
        % Compute derivatives of slack pressure wrt state variables
        [dpdx]     = GAS_dpdx(rho_ns,phi_0,phi_L,Funcs.SPF,Funcs.SFD,Funcs.NSFD);
    end
    
    % Call ddot (no conservation law is deleted)
    ddot = Inputs.MFd_CF;
    
    % Compute G vector for state variable derivative computation. This
    % vector does *not* include the slack-node line
    Gt1     = ddot;
    Gt2     = phi_0-phi_L;
    nsb_ind = [1:(Bus.Slack-1) (Bus.Slack+1):n]; % No slack bus
    nsl_ind = [1:(Bus.Ln_slk-1) (Bus.Ln_slk+1):m]; % No slack line
    G3p     = Gam.G3(nsl_ind,nsb_ind);
    G4p     = Gam.G4(nsl_ind,nsl_ind);
    G5p     = Gam.G5(nsl_ind,nsb_ind);
    Gt3p    = G3p*rho_ns - G4p*f(phi_L(nsl_ind),phi_0(nsl_ind),G5p*rho_ns);
    
    % Assemble
    Gt = [Gt1; Gt2; Gt3p];
    
    % Solve linear system [as a future reference, Mats.PC is an excellent
    % preconditioner which may be used to speed up inv(MM), if needed,
    % used approximately via dxdt_DAE = (Mats.PC*MM)\(Mats.PC*[Gt; 0]);]
    MM   = [Mats.Mx3  Mats.Mx4; 1 -dpdx];
    dxdt = MM\[Gt; 0];
    
    % The top variable is the slack-pressure derivative: set derivative output
    dxdt = dxdt(2:end);

end
end

% Local function: f
function vec = f(phi_L,phi_0,Gam5_rho)
    % Mass flux conservation
    vec = (phi_L+phi_0).*abs(phi_L+phi_0)./Gam5_rho;
end
