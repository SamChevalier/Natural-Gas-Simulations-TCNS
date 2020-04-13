function [h1,h2,Phi,Rho] = GAS_Build_Sig(gamma,phi_max,rho_nom,phi_MAX)
% GAS_BUILD_SIG: Build the sigmoid function and its derivative
%
% Input:
% 1) gamma     Speed of density drop-off
% 2) phi_max   Maximum flux flow of slack
% 3) rho_nom   Nominal density of slack. This isn't necessarily the stead
%              state value. It is the value which, when multiplied by S2,
%              gives the stead state density value
% 4) phi_MAX   The flux value at which the pressure drops by half
%
% Output:
% 1) h1        Sigmoid function derivative for S1
% 2) h2        Sigmoid function derivative for S2
% 3) Phi       Function for computing flux at the slack bus
% 4) Rho       Function for computing density at the slack bus
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms z phi_0

% Build the sigmoid functions
S1 = exp(z)./(1+exp(z));
S2 = exp(gamma*(phi_MAX - phi_0))./(1+exp(gamma*(phi_MAX - phi_0)));

% Build their derivatives
h1s = phi_max*(S1 - S1^2);
h2s = rho_nom*gamma*(S2^2 - S2);

% Flux and density functions
Phi_s = phi_max*S1;
Rho_s = rho_nom*S2;

% Turn into output functions
Phi = matlabFunction(Phi_s);
Rho = matlabFunction(Rho_s);
h1  = matlabFunction(h1s);
h2  = matlabFunction(h2s);

end
