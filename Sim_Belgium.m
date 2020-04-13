%% Simulate the Belgium System
%  26 Nodes; 27 Lines
clc; clear variables

% The values of C_ij^2 assume the following:
% flow     = 10^6 m^3/day
% pressure = bar (10^5 Pa)

% The flow and pressure units have been updated to the following:
% flow     = m^3/s
% pressure = Pa

%% A. Load System Data
ngc = Belgian_Sys_Data;  % "ngc" is a structure with elements:
                         % node, source, pipeline, compressor
% Bus and line counts
m = size(ngc.pipeline,1); % # of lines
n = size(ngc.node,1);     % # of buses

% Vector of nodal tpyes: 1 = const pres; 2 = const flux
vt                     = zeros(n,1);
vt(ngc.node(:,2) == 2) = 1;
vt(ngc.node(:,2) == 1) = 2;

% Source and load information
ind_vs      = find(vt ==1);
ind_vl      = find(vt ==2);
n_lds       = length(ind_vl);
source_pres = ngc.node(ind_vs,6);   % Fixed pressure at the sources
load_flux   = -ngc.node(ind_vl,3);  % Fixed flux at the loads

% Compressor ratio vector
cv                      = ones(n,1);
cv(ngc.compressor(:,1)) = ngc.compressor(:,2);

% Incidence matrix
E = GAS_Incidence(ngc);

% Build bus structure
Bus.Inc_M     = E;
Bus.Types     = vt;
Bus.Src_Press = source_pres;
Bus.Load_Injs = load_flux;
Bus.Comp_Rat  = cv;

%% B. Define Related Line and Friction Parameters

% Compute wave propagation
Line.Length   = ngc.pipeline(:,7);
Line.Diam     = ngc.pipeline(:,6);
Z             = 0.8;            % Compressibility
R             = 8.314462618;    % Ideal Gas Constant
T             = 281.15;         % Temperature
M             = 0.01900;        % Molecular Mass (NG)
Line.a        = sqrt(Z*R*T/M);  % Propagation speed
Line.a        = Line.a*ones(m,1);

% Compute Lambda - Friction Factor
ep            = 0.05e-3;                                % Absolute rugosity of pipe [m]
Line.Diam     = ngc.pipeline(:,6);                      % Line diameter [m]
Line.Lambda   = 1./((2*log10(3.7*Line.Diam/ep)).^2);    % Darcy Friction (unitless)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Alternatively, we may compute lambda via:          %
% % % Line.Lambda   = (C_ij^2)*D^5/(T*L*delta*Z*Cij^2) % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Predict the value of C_ij^2 - doesn't match exactly
C_ij2 = (Line.Diam)./(Line.Length.*Line.Lambda.*Line.a.^2);

%% C. Build Line Structure

% Optional initial guesses?
Bus.PG         = ngc.node(:,6);
Bus.DG         = -ngc.node(:,3);
Bus.DG(ind_vs) = ngc.source(:,2);
Line.FG        = ngc.pipeline(:,5);

% !!! Convert all pressure variables to densities !!! %
Bus.PG        = Bus.PG./(mean(Line.a)^2);             % Average a
Bus.Src_Press = Bus.Src_Press./(mean(Line.a)^2);      % Average a

% !!! Convert all flux parameters from m^3/s to kg/(s*m^2) !!! %
%               ~ use given density values ~
% Line.Area     = pi*(Line.Diam/2).^2;
% Bus.Load_Injs = Bus.Load_Injs.*Bus.PG(vt == 2)./(mean(Line.Area));
% Bus.DG        = Bus.DG.*Bus.PG./(mean(Line.Area));
% Line.FG       = Line.FG.*mean(Bus.PG)/mean(Line.Area);    % Just for initialization

%%%%%%%% Intialize the system via NR %%%%%%%%
NR.max_its      = 100;
NR.epsilon      = 0.0001;
[phi,rho,d_inj] = GAS_NR_Init(Bus,Line,NR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% D. Perform System Discretization & Re-Run Newton
DD             = 5e3;
[Bus_D,Line_D] = GAS_Sys_Discretization(Bus,Line,DD);

% Re-run NR to Initialize Discretized System
NR.max_its      = 100;
NR.epsilon      = 0.0001;
[phi,rho,d_inj] = GAS_NR_Init(Bus_D,Line_D,NR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Record Compressor Ratios
Comp_Rat = Bus_D.Comp_Rat;

% New m & n
m = size(Bus_D.Inc_M,1);
n = size(Bus_D.Inc_M,2);

% Set up simulation in puts (flux sources and loads)
Inputs.MFd_CP = zeros(n,1);  % Mass flux derivatives
Inputs.MFd_CF = zeros(n,1);  % Mass flux derivatives
Inputs.MFd_SF = zeros(n,1);  % Mass flux derivatives (one is removed)

%% 1) Test 1 - Load Increase at Node 16: Constant Pressure
Bus_D.Slack_Type = 'CP';
Bus_D.Slack      = 8;
Bus_D.CR_update  = 0;
Bus_D.Comp_Rat   = Comp_Rat;

% Define SS and ICs
phiL              = phi;
phi0              = phi;
Bus_D.SS          = [phi0; phiL; rho; d_inj];
IC                = Bus_D.SS;
Inputs.MFd_CP(16) = -0.125463;

% Simulation (a)
TD_Sim.tf = 1000;
TD_Sim.dt = 0.5;
[r1a]     = GAS_Dyn_Sim(IC,Bus_D,Line_D,Inputs,TD_Sim);

% Update ICs
Inputs.MFd_CP(16) = 0;
IC                = [r1a.phi0(end,:)'; r1a.phiL(end,:)'; r1a.rho(end,:)'; d_inj];

% Simulation (b)
TD_Sim.tf = 36*3600 - 1000;
TD_Sim.dt = 0.5;
[r1b]     = GAS_Dyn_Sim(IC,Bus_D,Line_D,Inputs,TD_Sim);

% Aggregate results
r1.rho  = [r1a.rho;  r1b.rho(2:end,:)];
r1.phi0 = [r1a.phi0; r1b.phi0(2:end,:)];
r1.phiL = [r1a.phiL; r1b.phiL(2:end,:)];
r1.t    = [r1a.t;    r1a.t(end) + r1b.t(2:end)];

%% 2) Test 2 - Load Increase at Node 16: Sigmoid Function
Bus_D.Slack_Type = 'SF';
Bus_D.Slack      = 8;
Bus_D.gamma      = 0.1;
Bus_D.phi_MAX    = d_inj(8)+100;
Bus_D.phi_max    = d_inj(8)+75;
Bus_D.CR_update  = 0;
Bus_D.Comp_Rat   = Comp_Rat;

% Define SS and ICs
phiL              = phi;
phi0              = phi;
Bus_D.SS          = [phi0; phiL; rho; d_inj];
IC                = Bus_D.SS;                 
Inputs.MFd_SF(16) = -0.125463;

% Simulation (a)
TD_Sim.tf = 1000;
TD_Sim.dt = 0.5;
[r2a]     = GAS_Dyn_Sim(IC,Bus_D,Line_D,Inputs,TD_Sim);

% Update ICs
Inputs.MFd_SF(16) = 0;
IC                = [r2a.phi0(end,:)'; r2a.phiL(end,:)'; r2a.rho(end,:)'; d_inj];

% Simulation (b)
TD_Sim.tf = 36*3600 - 1000;
TD_Sim.dt = 0.5;
[r2b]     = GAS_Dyn_Sim(IC,Bus_D,Line_D,Inputs,TD_Sim);

% Aggregate results
r2.rho  = [r2a.rho;  r2b.rho(2:end,:)];
r2.phi0 = [r2a.phi0; r2b.phi0(2:end,:)];
r2.phiL = [r2a.phiL; r2b.phiL(2:end,:)];
r2.z    = [r2a.z;    r2b.z(2:end,:)];
r2.t    = [r2a.t;    r2a.t(end) + r2b.t(2:end)];

%% 3) Test 3 - Load Increase at Node 16: Constant Flux
Bus_D.Slack_Type = 'CF';
Bus_D.Slack      = 8;
Bus_D.CR_update  = 0;
Bus_D.Comp_Rat   = Comp_Rat;

% Define SS and ICs
phiL              = phi;
phi0              = phi;
Bus_D.SS          = [phi0; phiL; rho; d_inj];
IC                = Bus_D.SS;
Inputs.MFd_CF(16) = -0.125463;

% Simulation (a)
TD_Sim.tf = 1000;
TD_Sim.dt = 0.5;
[r3a]     = GAS_Dyn_Sim(IC,Bus_D,Line_D,Inputs,TD_Sim);

% Update ICs
Inputs.MFd_CF(16) = 0;
IC                = [r3a.phi0(end,:)'; r3a.phiL(end,:)'; r3a.rho(end,:)'; d_inj];

% Simulation (b)
TD_Sim.tf = 36*3600 - 1000;
TD_Sim.dt = 0.5;
[r3b]     = GAS_Dyn_Sim(IC,Bus_D,Line_D,Inputs,TD_Sim);

% Aggregate results
r3.rho  = [r3a.rho;  r3b.rho(2:end,:)];
r3.phi0 = [r3a.phi0; r3b.phi0(2:end,:)];
r3.phiL = [r3a.phiL; r3b.phiL(2:end,:)];
r3.t    = [r3a.t;    r3a.t(end) + r3b.t(2:end)];

%% 4) Test 4 - Compressor Loss at Node 8: Constant Pressure
Bus_D.Slack_Type     = 'CP';
Bus_D.Slack          = 8;
Bus_D.Comp_Rat       = Comp_Rat;
Bus_D.CR_update      = 1;
Bus_D.CR_factors     = zeros(n,1);
Bus_D.CR_factors(17) = -0.03028;
Bus_D.CR_a           = 0.1;
Bus_D.CR_ts          = 100;

% Define SS and ICs
phiL     = phi;
phi0     = phi;
Bus_D.SS = [phi0; phiL; rho; d_inj];
IC       = Bus_D.SS;

% Simulation (a)
TD_Sim.tf = 200;
TD_Sim.dt = 0.5;
[r4a]     = GAS_Dyn_Sim(IC,Bus_D,Line_D,Inputs,TD_Sim);

% Update ICs
Bus_D.CR_update    = 0;
Bus_D.Comp_Rat(17) = 1.3028 - 0.03028;
IC                 = [r4a.phi0(end,:)'; r4a.phiL(end,:)'; r4a.rho(end,:)'; d_inj];

% Simulation (b)
TD_Sim.tf = 5*3600 - 200;
TD_Sim.dt = 0.5;
[r4b]     = GAS_Dyn_Sim(IC,Bus_D,Line_D,Inputs,TD_Sim);

% Aggregate results
r4.rho  = [r4a.rho;  r4b.rho(2:end,:)];
r4.phi0 = [r4a.phi0; r4b.phi0(2:end,:)];
r4.phiL = [r4a.phiL; r4b.phiL(2:end,:)];
r4.t    = [r4a.t;    r4a.t(end) + r4b.t(2:end)];

%% 5) Test 5 - Compressor Loss at Node 8: Sigmoid Function
Bus_D.Slack_Type     = 'SF';
Bus_D.Slack          = 8;
Bus_D.gamma          = 0.1;
Bus_D.phi_MAX        = d_inj(8)+100;
Bus_D.phi_max        = d_inj(8)+75;
Bus_D.Comp_Rat       = Comp_Rat;
Bus_D.CR_update      = 1;
Bus_D.CR_factors     = zeros(n,1);
Bus_D.CR_factors(17) = -0.03028;
Bus_D.CR_a           = 0.1;
Bus_D.CR_ts          = 100;

% Define SS and ICs
phiL     = phi;
phi0     = phi;
Bus_D.SS = [phi0; phiL; rho; d_inj];
IC       = Bus_D.SS;

% Simulation (a)
TD_Sim.tf = 200;
TD_Sim.dt = 0.5;
[r5a]     = GAS_Dyn_Sim(IC,Bus_D,Line_D,Inputs,TD_Sim);

% Update ICs
Bus_D.CR_update    = 0;
Bus_D.Comp_Rat(17) = 1.3028 - 0.03028;
IC                 = [r5a.phi0(end,:)'; r5a.phiL(end,:)'; r5a.rho(end,:)'; d_inj];

% Simulation (b)
TD_Sim.tf = 5*3600 - 200;
TD_Sim.dt = 0.5;
[r5b]     = GAS_Dyn_Sim(IC,Bus_D,Line_D,Inputs,TD_Sim);

% Aggregate results
r5.rho  = [r5a.rho;  r5b.rho(2:end,:)];
r5.phi0 = [r5a.phi0; r5b.phi0(2:end,:)];
r5.phiL = [r5a.phiL; r5b.phiL(2:end,:)];
r5.t    = [r5a.t;    r5a.t(end) + r5b.t(2:end)];

%% 6) Test 6 - Compressor Loss at Node 8: Constant Flux
Bus_D.Slack_Type     = 'CF';
Bus_D.Slack          = 8;
Bus_D.CR_update      = 0;
Bus_D.Comp_Rat       = Comp_Rat;
Bus_D.CR_update      = 1;
Bus_D.CR_factors     = zeros(n,1);
Bus_D.CR_factors(17) = -0.03028;
Bus_D.CR_a           = 0.1;
Bus_D.CR_ts          = 100;

% Define SS and ICs
phiL     = phi;
phi0     = phi;
Bus_D.SS = [phi0; phiL; rho; d_inj];
IC       = Bus_D.SS;

% Simulation (a)
TD_Sim.tf = 200;
TD_Sim.dt = 0.5;
[r6a]     = GAS_Dyn_Sim(IC,Bus_D,Line_D,Inputs,TD_Sim);

% Update ICs
Bus_D.CR_update    = 0;
Bus_D.Comp_Rat(17) = 1.3028 - 0.03028;
IC                 = [r6a.phi0(end,:)'; r6a.phiL(end,:)'; r6a.rho(end,:)'; d_inj];

% Simulation (b)
TD_Sim.tf = 5*3600 - 200;
TD_Sim.dt = 0.5;
[r6b]     = GAS_Dyn_Sim(IC,Bus_D,Line_D,Inputs,TD_Sim);

% Aggregate results
r6.rho  = [r6a.rho;  r6b.rho(2:end,:)];
r6.phi0 = [r6a.phi0; r6b.phi0(2:end,:)];
r6.phiL = [r6a.phiL; r6b.phiL(2:end,:)];
r6.t    = [r6a.t;    r6a.t(end) + r6b.t(2:end)];

%% Colors
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250];
c4 = [0.4940, 0.1840, 0.5560];
c5 = [0.4660, 0.6740, 0.1880];
c6 = [0.3010, 0.7450, 0.9330];
c7 = [0.6350, 0.0780, 0.1840];

%% Plot Test 1: Flux Zoom
close all
stp = 16000;
plot(r1.t(1:25:stp)/3600,r1.phiL((1:25:stp),[1:37 39:70 72:end]),'Linewidth',0.1)
hold on
a = plot(r1.t(1:25:stp)/3600,r1.phi0((1:25:stp),8 ),'-.','Linewidth',3,'color',c1);
b = plot(r1.t(1:25:stp)/3600,r1.phi0((1:25:stp),38),':','Linewidth',3,'color',c2);
c = plot(r1.t(1:25:stp)/3600,r1.phiL((1:25:stp),71),'-','Linewidth',3,'color',c6);
legend([a, b, c],{'$\phi_5$','$\phi_8$','$\phi_{16}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'$\rm{Gas \; Flux \; \phi \;(\frac{kg}{s\cdot m^2})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
xlim([0 1])
ylim([-50 300])
xticks(0:0.2:1)
grid on

% Size
set(gcf,'Units','inches','Position',[0 0 9 4])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% Plot Test 1: Flux Full
close all
ind = 100;
plot(r1.t(1:ind:end)/3600,r1.phiL((1:ind:end),[1:37 39:70 72:end]),'Linewidth',0.1)
hold on
a = plot(r1.t(1:ind:end)/3600,r1.phi0((1:ind:end),8 ),'-.','Linewidth',3,'color',c1);
b = plot(r1.t(1:ind:end)/3600,r1.phi0((1:ind:end),38),':','Linewidth',3,'color',c2);
c = plot(r1.t(1:ind:end)/3600,r1.phiL((1:ind:end),71),'-','Linewidth',3,'color',c6);
legend([a, b, c],{'$\phi_5$','$\phi_8$','$\phi_{16}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'$\rm{Gas \; Flux \; \phi \;(\frac{kg}{s\cdot m^2})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
xlim([0 36])
ylim([-50 420])
xticks(0:5:35)
grid on

% Size
set(gcf,'Units','inches','Position',[0 0 9 4])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% Plot Test 1: Density Full
close all
ind = 100;
plot(r1.t(1:ind:end)/3600,r1.rho((1:ind:end),[1:37 39:70 72:end]),'Linewidth',0.1)
hold on
a = plot(r1.t(1:ind:end)/3600,r1.rho((1:ind:end),5 ),'-.','Linewidth',3,'color',c1);
b = plot(r1.t(1:ind:end)/3600,r1.rho((1:ind:end),8),':','Linewidth',3,'color',c2);
c = plot(r1.t(1:ind:end)/3600,r1.rho((1:ind:end),16),'-','Linewidth',3,'color',c6);
legend([a, b, c],{'$\rho_5$','$\rho_8$','$\rho_{16}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'$\rm{Gas \; Density \; \rho \;(\frac{kg}{m^3})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
xlim([0 36])
ylim([35 77])
xticks(0:5:35)
grid on

% Size
set(gcf,'Units','inches','Position',[0 0 9 4])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% (**Paper Plot**) Plot Test 1: Flux and Density Full
close all
subplot(1,2,1)
ind = 1000;
plot(r1.t(1:ind:end)/3600,r1.phiL((1:ind:end),[1:37 39:70 72:end]),'Linewidth',0.1)
hold on
a = plot(r1.t(1:ind:end)/3600,r1.phi0((1:ind:end),8 ),'-.','Linewidth',3,'color',c1);
b = plot(r1.t(1:ind:end)/3600,r1.phi0((1:ind:end),38),':','Linewidth',3,'color',c2);
c = plot(r1.t(1:ind:end)/3600,r1.phiL((1:ind:end),71),'-','Linewidth',3,'color',c6);
legend([a, b, c],{'$\phi_5$','$\phi_8$','$\phi_{16}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'$\rm{Gas \; Flux \; \phi \;(\frac{kg}{s\cdot m^2})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
text(0.95,400,'$({\bf a})$','Interpreter','latex','FontSize',13)
xlim([0 36])
ylim([-50 420])
xticks(0:10:30)
grid on

subplot(1,2,2)
ind = 100;
plot(r1.t(1:ind:end)/3600,r1.rho((1:ind:end),[1:37 39:70 72:end]),'Linewidth',0.1)
hold on
a = plot(r1.t(1:ind:end)/3600,r1.rho((1:ind:end),5 ),'-.','Linewidth',3,'color',c1);
b = plot(r1.t(1:ind:end)/3600,r1.rho((1:ind:end),8),':','Linewidth',3,'color',c2);
c = plot(r1.t(1:ind:end)/3600,r1.rho((1:ind:end),16),'-','Linewidth',3,'color',c6);
legend([a, b, c],{'$\rho_5$','$\rho_8$','$\rho_{16}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'$\rm{Gas \; Density \; \rho \;(\frac{kg}{m^3})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
text(0.95,77.8,'$({\bf b})$','Interpreter','latex','FontSize',13)
xlim([0 36])
ylim([36 80])
xticks(0:10:30)
yticks(40:10:80)
grid on

% Size
set(gcf,'Units','inches','Position',[0 0 9 4])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% (**Paper Plot**) Plot Test 1: Flux and Density Full - Top/Bottom
close all
subplot(2,1,1)
ind = 1000;
plt = plot(r1.t(1:ind:end)/3600,r1.phiL((1:ind:end),[1:37 39:70 72:end]),'Linewidth',0.1);
for ii = 1:length(plt)
    clr = plt(ii).Color;
    plt(ii).Color = [clr 0.2];
end
hold on
a = plot(r1.t(1:ind:end)/3600,r1.phi0((1:ind:end),8 ),'-.','Linewidth',3,'color',c1);
b = plot(r1.t(1:ind:end)/3600,r1.phi0((1:ind:end),38),':','Linewidth',3,'color',c2);
c = plot(r1.t(1:ind:end)/3600,r1.phiL((1:ind:end),71),'-','Linewidth',3,'color',c6);

numcolumns = 3;
str        = {'$\phi_5$','$\phi_8$','$\phi_{16}$'};
plots      = [a, b, c];
[legend_h,object_h,plot_h,text_strings] = columnlegend(numcolumns,str,plots);

ylabel({'$\rm{Gas \; Flux \; \phi \;(\frac{kg}{s\cdot m^2})}$'},'Interpreter','latex','FontSize',14)
%xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
text(0.4,440,'$({\bf a})$','Interpreter','latex','FontSize',13)
xlim([0 36])
ylim([-50 500])
xticks(0:10:30)
yticks([0 200 400])
grid on

subplot(2,1,2)
plt = plot(r1.t(1:ind:end)/3600,r1.rho((1:ind:end),[1:37 39:70 72:end]),'Linewidth',0.1);
for ii = 1:length(plt)
    clr = plt(ii).Color;
    plt(ii).Color = [clr 0.2];
end
hold on
a = plot(r1.t(1:ind:end)/3600,r1.rho((1:ind:end),5 ),'-.','Linewidth',3,'color',c1);
b = plot(r1.t(1:ind:end)/3600,r1.rho((1:ind:end),8),':','Linewidth',3,'color',c2);
c = plot(r1.t(1:ind:end)/3600,r1.rho((1:ind:end),16),'-','Linewidth',3,'color',c6);

numcolumns = 3;
str        = {'$\rho_5$','$\rho_8$','$\rho_{16}$'};
plots      = [a, b, c];
[legend_h,object_h,plot_h,text_strings] = columnlegend(numcolumns,str,plots);

ylabel({'$\rm{Gas \; Density \; \rho \;(\frac{kg}{m^3})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
text(0.4,80,'$({\bf b})$','Interpreter','latex','FontSize',13)
xlim([0 36])
ylim([36 85])
xticks(0:10:30)
yticks(40:10:80)
grid on

% Size
set(gcf,'Units','inches','Position',[0 0 9 5])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% (**Paper Plot**) Plot Test 2: Flux and Density Full
close all
subplot(1,2,1)
ind = 1000;
plot(r2.t(1:ind:end)/3600,r2.phiL((1:ind:end),[1:37 39:70 72:end]),'Linewidth',0.1)
hold on
a = plot(r2.t(1:ind:end)/3600,r2.phi0((1:ind:end),8 ),'-.','Linewidth',3,'color',c1);
b = plot(r2.t(1:ind:end)/3600,r2.phi0((1:ind:end),38),':','Linewidth',3,'color',c2);
c = plot(r2.t(1:ind:end)/3600,r2.phiL((1:ind:end),71),'-','Linewidth',3,'color',c6);
legend([a, b, c],{'$\phi_5$','$\phi_8$','$\phi_{16}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'$\rm{Gas \; Flux \; \phi \;(\frac{kg}{s\cdot m^2})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
text(0.95,400,'$({\bf a})$','Interpreter','latex','FontSize',13)
xlim([0 36])
ylim([-50 420])
xticks(0:10:30)
grid on

subplot(1,2,2)
ind = 100;
plot(r2.t(1:ind:end)/3600,r2.rho((1:ind:end),[1:37 39:70 72:end]),'Linewidth',0.1)
hold on
a = plot(r2.t(1:ind:end)/3600,r2.rho((1:ind:end),5 ),'-.','Linewidth',3,'color',c1);
b = plot(r2.t(1:ind:end)/3600,r2.rho((1:ind:end),8),':','Linewidth',3,'color',c2);
c = plot(r2.t(1:ind:end)/3600,r2.rho((1:ind:end),16),'-','Linewidth',3,'color',c6);
legend([a, b, c],{'$\rho_5$','$\rho_8$','$\rho_{16}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'$\rm{Gas \; Density \; \rho \;(\frac{kg}{m^3})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
text(0.95,77.8,'$({\bf b})$','Interpreter','latex','FontSize',13)
xlim([0 36])
ylim([28 80])
xticks(0:10:30)
grid on

% Size
set(gcf,'Units','inches','Position',[0 0 9 4])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% (**Paper Plot**) Plot Test 2: Flux and Density Full - Top/Bottom
close all
subplot(2,1,1)
ind = 1000;
plt = plot(r2.t(1:ind:end)/3600,r2.phiL((1:ind:end),[1:37 39:70 72:end]),'Linewidth',0.1);
for ii = 1:length(plt)
    clr = plt(ii).Color;
    plt(ii).Color = [clr 0.2];
end
hold on
a = plot(r2.t(1:ind:end)/3600,r2.phi0((1:ind:end),8 ),'-.','Linewidth',3,'color',c1);
b = plot(r2.t(1:ind:end)/3600,r2.phi0((1:ind:end),38),':','Linewidth',3,'color',c2);
c = plot(r2.t(1:ind:end)/3600,r2.phiL((1:ind:end),71),'-','Linewidth',3,'color',c6);

numcolumns = 3;
str        = {'$\phi_5$','$\phi_8$','$\phi_{16}$'};
plots      = [a, b, c];
[legend_h,object_h,plot_h,text_strings] = columnlegend(numcolumns,str,plots);

ylabel({'$\rm{Gas \; Flux \; \phi \;(\frac{kg}{s\cdot m^2})}$'},'Interpreter','latex','FontSize',14)
%xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
text(0.4,445,'$({\bf a})$','Interpreter','latex','FontSize',13)
xlim([0 36])
ylim([-50 500])
xticks(0:10:30)
yticks([0 200 400])
grid on

subplot(2,1,2)
plt = plot(r2.t(1:ind:end)/3600,r2.rho((1:ind:end),[1:37 39:70 72:end]),'Linewidth',0.1);
for ii = 1:length(plt)
    clr = plt(ii).Color;
    plt(ii).Color = [clr 0.2];
end
hold on
a = plot(r2.t(1:ind:end)/3600,r2.rho((1:ind:end),5 ),'-.','Linewidth',3,'color',c1);
b = plot(r2.t(1:ind:end)/3600,r2.rho((1:ind:end),8),':','Linewidth',3,'color',c2);
c = plot(r2.t(1:ind:end)/3600,r2.rho((1:ind:end),16),'-','Linewidth',3,'color',c6);

numcolumns = 3;
str        = {'$\rho_5$','$\rho_8$','$\rho_{16}$'};
plots      = [a, b, c];
[legend_h,object_h,plot_h,text_strings] = columnlegend(numcolumns,str,plots);

ylabel({'$\rm{Gas \; Density \; \rho \;(\frac{kg}{m^3})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
text(0.4,80,'$({\bf b})$','Interpreter','latex','FontSize',13)
xlim([0 36])
ylim([28 85])
xticks(0:10:30)
yticks(40:10:80)
grid on

% Size
set(gcf,'Units','inches','Position',[0 0 9 5])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% Plot Test 2: z
close all
plot(r2.t(1:25:end)/3600,r2.z(1:25:end),'Linewidth',1)
ylabel({'${\rm State} \; {\rm Variable}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
legend({'$z$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
set(gca,'FontName','Times','FontSize',14)
xlim([0 36])
ylim([0 11])
xticks(0:10:30)
grid on

% Size
set(gcf,'Units','inches','Position',[0 0 9 2])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% (**Paper Plot**) Plot Test 3: Flux and Density Full
close all
subplot(1,2,1)
ind = 1000;
plot(r3.t(1:ind:end)/3600,r3.phiL((1:ind:end),[1:37 39:70 72:end]),'Linewidth',0.1)
hold on
a = plot(r3.t(1:ind:end)/3600,r3.phi0((1:ind:end),8 ),'-.','Linewidth',3,'color',c1);
b = plot(r3.t(1:ind:end)/3600,r3.phi0((1:ind:end),38),':','Linewidth',3,'color',c2);
c = plot(r3.t(1:ind:end)/3600,r3.phiL((1:ind:end),71),'-','Linewidth',3,'color',c6);
legend([a, b, c],{'$\phi_5$','$\phi_8$','$\phi_{16}$'},'Interpreter','latex','box','off','Location','NE','FontSize',12)
ylabel({'$\rm{Gas \; Flux \; \phi \;(\frac{kg}{s\cdot m^2})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
text(0.95,315,'$({\bf a})$','Interpreter','latex','FontSize',13)
xlim([0 36])
ylim([-50 335])
xticks(0:10:30)
grid on

subplot(1,2,2)
plot(r3.t(1:ind:end)/3600,r3.rho((1:ind:end),[1:37 39:70 72:end]),'Linewidth',0.1)
hold on
a = plot(r3.t(1:ind:end)/3600,r3.rho((1:ind:end),5 ),'-.','Linewidth',3,'color',c1);
b = plot(r3.t(1:ind:end)/3600,r3.rho((1:ind:end),8),':','Linewidth',3,'color',c2);
c = plot(r3.t(1:ind:end)/3600,r3.rho((1:ind:end),16),'-','Linewidth',3,'color',c6);
legend([a, b, c],{'$\rho_5$','$\rho_8$','$\rho_{16}$'},'Interpreter','latex','box','off','Location','NE','FontSize',12)
ylabel({'$\rm{Gas \; Density \; \rho \;(\frac{kg}{m^3})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
text(0.95,78.5,'$({\bf b})$','Interpreter','latex','FontSize',13)
xlim([0 36])
ylim([0 83])
xticks(0:10:30)
grid on

% Size
set(gcf,'Units','inches','Position',[0 0 9 4])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% (**Paper Plot**) Plot Test 3: Flux and Density Full - Top/Bottom
close all
subplot(2,1,1)
ind = 1000;
plt = plot(r3.t(1:ind:end)/3600,r3.phiL((1:ind:end),[1:37 39:70 72:end]),'Linewidth',0.1);
for ii = 1:length(plt)
    clr = plt(ii).Color;
    plt(ii).Color = [clr 0.2];
end
hold on
a = plot(r3.t(1:ind:end)/3600,r3.phi0((1:ind:end),8 ),'-.','Linewidth',3,'color',c1);
b = plot(r3.t(1:ind:end)/3600,r3.phi0((1:ind:end),38),':','Linewidth',3,'color',c2);
c = plot(r3.t(1:ind:end)/3600,r3.phiL((1:ind:end),71),'-','Linewidth',3,'color',c6);

numcolumns = 3;
str        = {'$\phi_5$','$\phi_8$','$\phi_{16}$'};
plots      = [a, b, c];
[legend_h,object_h,plot_h,text_strings] = columnlegend(numcolumns,str,plots);

ylabel({'$\rm{Gas \; Flux \; \phi \;(\frac{kg}{s\cdot m^2})}$'},'Interpreter','latex','FontSize',14)
%xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
text(0.4,360,'$({\bf a})$','Interpreter','latex','FontSize',13)
xlim([0 36])
ylim([-50 400])
xticks(0:10:30)
grid on

subplot(2,1,2)
plt = plot(r3.t(1:ind:end)/3600,r3.rho((1:ind:end),[1:37 39:70 72:end]),'Linewidth',0.1);
for ii = 1:length(plt)
    clr = plt(ii).Color;
    plt(ii).Color = [clr 0.2];
end
hold on
a = plot(r3.t(1:ind:end)/3600,r3.rho((1:ind:end),5 ),'-.','Linewidth',3,'color',c1);
b = plot(r3.t(1:ind:end)/3600,r3.rho((1:ind:end),8),':','Linewidth',3,'color',c2);
c = plot(r3.t(1:ind:end)/3600,r3.rho((1:ind:end),16),'-','Linewidth',3,'color',c6);

numcolumns = 3;
str        = {'$\rho_5$','$\rho_8$','$\rho_{16}$'};
plots      = [a, b, c];
[legend_h,object_h,plot_h,text_strings] = columnlegend(numcolumns,str,plots);

ylabel({'$\rm{Gas \; Density \; \rho \;(\frac{kg}{m^3})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
text(0.4,82,'$({\bf b})$','Interpreter','latex','FontSize',13)
xlim([0 36])
ylim([0 90])
xticks(0:10:30)
grid on

% Size
set(gcf,'Units','inches','Position',[0 0 9 5])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% (**Paper Plot**) Plot Test 4 and 6: Flux Reverse Flows
close all
ind = 10;
subplot(1,2,1)
plot(r4.t(1:ind:end)/3600,r4.phi0((1:ind:end),79:end),'Linewidth',0.65)
xlim([0 0.6])
ylabel({'$\rm{Gas \; Flux \; \phi \;(\frac{kg}{s\cdot m^2})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
grid on
text(0.015,33,'$({\bf a})$','Interpreter','latex','FontSize',13)

subplot(1,2,2)
plot(r6.t(1:ind:end)/3600,r6.phi0((1:ind:end),79:end),'Linewidth',0.65)
xlim([0 0.6])
ylabel({'$\rm{Gas \; Flux \; \phi \;(\frac{kg}{s\cdot m^2})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
grid on
text(0.015,33,'$({\bf b})$','Interpreter','latex','FontSize',13)

% Size
set(gcf,'Units','inches','Position',[0 0 9 3])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% (**Paper Plot**) Plot Test 4 and 6: Flux Reverse Flows - Top/Bottom
close all
ind = 10;
subplot(2,1,1)
plot(r4.t(1:ind:end)/3600,r4.phi0((1:ind:end),79:end),'Linewidth',0.65)
xlim([0 0.6])
ylabel({'$\rm{Gas \; Flux \; \phi \;(\frac{kg}{s\cdot m^2})}$'},'Interpreter','latex','FontSize',14)
%xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
grid on
text(0.006,33,'$({\bf a})$','Interpreter','latex','FontSize',13)

subplot(2,1,2)
plot(r6.t(1:ind:end)/3600,r6.phi0((1:ind:end),79:end),'Linewidth',0.65)
xlim([0 0.6])
ylabel({'$\rm{Gas \; Flux \; \phi \;(\frac{kg}{s\cdot m^2})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
grid on
text(0.006,33,'$({\bf b})$','Interpreter','latex','FontSize',13)

% Size
set(gcf,'Units','inches','Position',[0 0 9 4])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% (**Paper Plot**) Plot Test 4 and 6: Nodal Pressure (zoom in to load nodes 18, 19, 20)
close all
ind = 10;
subplot(1,2,1)
plot(r4.t(1:ind:end)/3600,r4.rho((1:ind:end),:),'Linewidth',0.65)
xlim([0 5])
ylim([74.1 76])
ylabel({'$\rm{Gas \; Density \; \rho \;(\frac{kg}{m^3})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
grid on
text(4.5,75.87,'$({\bf a})$','Interpreter','latex','FontSize',13)

subplot(1,2,2)
plot(r6.t(1:ind:end)/3600,r6.rho((1:ind:end),:),'Linewidth',0.65)
xlim([0 5])
ylim([74.1 76])
ylabel({'$\rm{Gas \; Density \; \rho \;(\frac{kg}{m^3})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
grid on
text(4.5,75.87,'$({\bf b})$','Interpreter','latex','FontSize',13)

% Size
set(gcf,'Units','inches','Position',[0 0 9 3])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% (**Paper Plot**) Plot Test 4 and 6: Nodal Pressure (zoom in to load nodes 18, 19, 20) - Top/Bottom
close all
ind = 10;
subplot(2,1,1)
plot(r4.t(1:ind:end)/3600,r4.rho((1:ind:end),:),'Linewidth',0.65)
xlim([0 3])
ylim([74 76])
ylabel({'$\rm{Gas \; Density \; \rho \;(\frac{kg}{m^3})}$'},'Interpreter','latex','FontSize',14)
%xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
grid on
text(2.85,75.72,'$({\bf a})$','Interpreter','latex','FontSize',13)

subplot(2,1,2)
plot(r6.t(1:ind:end)/3600,r6.rho((1:ind:end),:),'Linewidth',0.65)
xlim([0 3])
ylim([74 76])
ylabel({'$\rm{Gas \; Density \; \rho \;(\frac{kg}{m^3})}$'},'Interpreter','latex','FontSize',14)
xlabel({'$\rm{Time \; (hours)}$'},'Interpreter','latex','FontSize',14)
set(gca,'FontName','Times','FontSize',14)
grid on
text(2.85,75.72,'$({\bf b})$','Interpreter','latex','FontSize',13)

% Size
set(gcf,'Units','inches','Position',[0 0 9 4.5])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"
