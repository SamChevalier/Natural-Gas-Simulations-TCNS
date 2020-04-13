function [VBus,VLine] = GAS_Sys_Discretization(Bus,Line,DD)
% GAS_Sys_Discretization:  Perform discretization of all lines and
%                          a set of virtual nodes.
%
% Input:
% 1) Bus       Node structure with the following elements
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
% 2) Line      Line structure with the following elements
%                   Length     => Directed incidence matrix
%                   Lambda     => 1=source, 2=load
%                   Diam       => Load (2) injections 
%                   a          => Propagation/speed factor
% 3) DD         Discretization distance
%
% Output:
% 1) p_s       Numerical value of the source pressure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bus and line counts
E = Bus.Inc_M;
m = size(E,1); % # of lines
n = size(E,2); % # of buses

% Initialize Virtual Bus Structures
VBus.Inc_M     = [];
VBus.Types     = Bus.Types;
VBus.Src_Press = Bus.Src_Press;  % Careful: these are listed in order
VBus.Load_Injs = Bus.Load_Injs;  % Careful: these are listed in order
VBus.Comp_Rat  = Bus.Comp_Rat;

% Initialize Virtual Line Structures
VLine.Length   = [];
VLine.Lambda   = [];
VLine.Diam     = [];
VLine.a        = [];

% Get node list
node_start = find(0.5*(E+abs(E))'==1) - (0:n:(n*(m-1)))';
node_end   = find(0.5*(abs(E)-E)'==1) - (0:n:(n*(m-1)))';

% Line list
V_lines = [];
V_nodes = 1:n;

% List of start and end nodes, correlated with lines (initially empty)
All_start_nodes = [];
All_end_nodes   = [];

% Initially, loop over all lines and segment appropriately
for ii = 1:m
    % If the division has a remainder, add an extra segment if the
    % remainder is 0.5 or larger. Otherwise, round down.
    num_Vlines = round(Line.Length(ii)/DD);
    if num_Vlines == 0
        num_Vlines = 1;
    end
    
    % How many new nodes? One less than the number of segments
    num_Vnodes = num_Vlines - 1;
    
    % What are the indices of the nodes these lines are tied to?
    new_Vnodes          = V_nodes(end) + (1:num_Vnodes);
    V_nodes(new_Vnodes) = new_Vnodes;
    
    % Define line indices
    if ii == 1
        line_inds          = 1:num_Vlines;
        V_lines(line_inds) = line_inds;
    else
        line_inds          = line_inds(end) + (1:num_Vlines);
        V_lines(line_inds) = line_inds;
    end
    % What are the characteristics of these new lines?
    VLine.Length(line_inds) = ones(num_Vlines,1)*Line.Length(ii)/num_Vlines;
    VLine.Lambda(line_inds) = ones(num_Vlines,1)*Line.Lambda(ii);
    VLine.Diam(line_inds)   = ones(num_Vlines,1)*Line.Diam(ii);
    VLine.a(line_inds)      = ones(num_Vlines,1)*Line.a(ii);
    
    % Update VBus structure
    ns                            = length(VBus.Src_Press);
    VBus.Types(new_Vnodes)        = 2; 
    VBus.Load_Injs(new_Vnodes-ns) = 0; % Subtract out the number of source
                                       % indices, which will necessarily be
                                       % behind the new virtual nodes
    VBus.Comp_Rat(new_Vnodes)     = 1;
    
    % Build list of nodes
    nodes                   = [node_start(ii) new_Vnodes node_end(ii)];  
    [start_nodes,end_nodes] = node_list(nodes);
    All_start_nodes         = [All_start_nodes start_nodes];
    All_end_nodes           = [All_end_nodes   end_nodes];
end

% Transpose
V_nodes         = V_nodes';
All_start_nodes = All_start_nodes';
All_end_nodes   = All_end_nodes';
VLine.Length    = VLine.Length';
VLine.Lambda    = VLine.Lambda';
VLine.Diam      = VLine.Diam';
VLine.a         = VLine.a';

% Construct a new Incidence Matrix
ngc.node     = V_nodes;
ngc.pipeline = [All_start_nodes All_end_nodes];
VBus.Inc_M   = GAS_Incidence(ngc);

end

% Local function: convert nodes
function [start_nodes,end_nodes] = node_list(nodes)
    % Convert a list of nodes into a list of start and end nodes
    nv          = length(nodes);
    start_nodes = nodes(1:nv-1);
    end_nodes   = nodes(2:nv);
end