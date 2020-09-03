%% --- Flow within an idealised zebrafish ISV network
%  Lowell Taylor Edgar
%  Usher Institute, University of Edinburgh
%  2020
%--------------------------------------------------------------------------
clc
clear all
close all

%% --- Input parameters
load('Ctrl_fish_240_stats.mat')

% Number of aISVs
NaISVs = 15;

% Number of vISVs
NvISVs = NaISVs;

% Dynamic viscosity of blood (kg/um/s)
mu = 3.5*1e-3*1e-6; 

% Dynamic viscosity of plasma/water (kg/um/s) (for GATA mutants)
%mu = 1.0*1e-3*1e-6;

% Define pressure boundary conditions (kg/s^2/um, 1e-6 Pa)
% Dorsal aorta inlet pressure
PDA1 = 201.3168 *1e-6;

% Length over which the pressure of the aorta is dissipated (um)
L_aorta_dissipate = 2000;

% Dorsal aorta outlet pressure (defined later with length of aorta)
%PDA2 = PDA1*0.5;
% Pectoral cardinal vein pressure
PPCV = 0;

% Numer of time steps
%TS = 8;

% Diameter of dorsal aorta (um)
%DDA = [18.2532   24.2864   25.0368   24.3244   23.1058   23.1464   23.9432   25.4234];

assert(length(DDA) == TS, "Not enough time steps in DA diameter data")

% Diameter of aISVs (um)
%DaISV = [6.8123    6.0493    5.8603    6.2328    7.0138    8.5477    7.1855    6.9223];
%DaISV = DaISV(1)*ones(1, TS);

assert(length(DaISV) == TS, "Not enough time steps in aISV diameter data")

% Diameter of vISVs (um)
%DvISV = [7.4278    5.7623    5.6581    6.4274    6.7766    7.8029    8.0949    7.5754];
%DvISV = sum([2*DvISV(1)/3 + DvISV(TS)/3])*ones(1, TS);

assert(length(DvISV) == TS, "Not enough time steps in vISV diameter data")

% Diameter of DLAV (um)
%DDLAV = [0.0000    0.0000    0.0000    5.4870    7.4895    8.1460    7.8604    8.0289];

assert(length(DDLAV) == TS, "Not enough time steps in DLAV diameter data")

% Flag to display network with numbering prior to solving for flow
plot_network_with_numbering = true;

%--------------------------------------------------------------------------
%% --- Create the network
nodes = [];
edges = [];
vess_type = [];

aorta_BC_nodes = [];
vISV_BC_nodes = [];

% --- Create the dorsal aorta
node1 = [-lsep/2 0];
node2 = [0 0];

nodes = [nodes; node1; node2];

edge1 = [1 2 norm(node2 - node1)];

edges = [edges; edge1];
vess_type = [vess_type; "aorta"];

aorta_BC_nodes = [aorta_BC_nodes; 1];

for i = 2:NaISVs
    node1 = length(nodes);
    old_node = nodes(node1,:);
    new_node = [old_node(1)+lsep 0];

    nodes = [nodes; new_node];
    node2 = length(nodes);
    
    new_edge = [node1 node2 norm(new_node - old_node)];
    edges = [edges; new_edge];
    vess_type = [vess_type; "aorta"];
end

node1 = length(nodes);
old_node = nodes(node1,:);

if (NvISVs >= NaISVs)
    last_node = [NvISVs*lsep 0];
else
    last_node = [NaISVs*lsep 0];
end

nodes = [nodes; last_node];
node2 = length(nodes);

last_edge = [node1 node2 norm(last_node - old_node)];
edges = [edges; last_edge];
vess_type = [vess_type; "aorta"];

end_aorta_node = length(nodes);
aorta_BC_nodes = [aorta_BC_nodes; end_aorta_node];

aorta_length = sum(edges(:,3));

% Define the outlet aortic pressure
PDA2 = PDA1*(1 - aorta_length/L_aorta_dissipate);
%PDA2 = PDA1;

% --- Create the aISVs
aISV_end_nodes = [];

for i = 1:NaISVs
    xpt = (i-1)*lsep;
    
    node1 = i+1;
    old_node = nodes(node1,:);
    
    for k = 1:3
        ypt = (k/3)*LISV;
        
        new_node = [xpt ypt];
        nodes = [nodes; new_node];
        node2 = length(nodes);
        
        new_edge = [node1 node2 norm(new_node - old_node)];
        edges = [edges; new_edge];
        vess_type = [vess_type; "aISV"];
        
        node1 = node2;
        old_node = nodes(node1,:);
    end
    
    aISV_end_nodes = [aISV_end_nodes; node2];
end

% --- Create the vISVs
vISV_end_nodes = [];

for i = 1:2:2*NvISVs
    xpt = i*lsep/2;
    
    old_node = [xpt -1];
    nodes = [nodes; old_node];
    node1 = length(nodes);
    
    vISV_BC_nodes = [vISV_BC_nodes; node1];
    
    for k = 1:3
        ypt = (k/3)*LISV;
        
        new_node = [xpt ypt];
        nodes = [nodes; new_node];
        node2 = length(nodes);
        
        new_edge = [node1 node2 norm(new_node - old_node)];
        edges = [edges; new_edge];
        vess_type = [vess_type; "vISV"];
        
        node1 = node2;
        old_node = nodes(node1,:);
    end
    
    vISV_end_nodes = [vISV_end_nodes; node1];
end

% --- Create the DLAV
% Create segment associated with aISV
for i = 1:NaISVs
    node1 = aISV_end_nodes(i);
    old_node = nodes(node1);
    
    if (i <= NvISVs)
        node2 = vISV_end_nodes(i);
    else
        if ((i+1) < NaISVs)
            node2 = aISV_end_nodes(i+1);
        end
    end
    
    new_node = nodes(node2);

    new_edge = [node1 node2 norm(new_node - old_node)];
    edges = [edges; new_edge];
    vess_type = [vess_type; "dlav"];
end

% Create segment associated with vISV
for i = 1:NvISVs-1
    node1 = vISV_end_nodes(i);
    old_node = nodes(node1,:);
    
    if ((i+1) <= NaISVs)
        node2 = aISV_end_nodes(i+1);
    else
        node2 = vISV_end_nodes(i+1);
    end
    
    new_node = nodes(node2,:);

    new_edge = [node1 node2 norm(new_node - old_node)];
    edges = [edges; new_edge];
    vess_type = [vess_type; "dlav"];
end       

num_nodes = length(nodes);
num_edges = length(edges);

% --- Create the table for diameter over time 
D = zeros(num_edges, TS);

for i = 1:num_edges
    if (vess_type(i) == "aorta")
        D(i,:) = DDA;
    end
    
    if (vess_type(i) == "aISV")
        D(i,:) = DaISV;
    end
    
    if (vess_type(i) == "vISV")
        D(i,:) = DvISV;
    end
    
    if (vess_type(i) == "dlav")
        D(i,:) = DDLAV;
    end
end
 
for v = 1:num_edges
    for t = 1:TS
        if (D(v,t) == 0)
            D(v,t) = 1e-10;
        end
    end
end

%--------------------------------------------------------------------------
%% --- Display the network configuration and numbering before running flow solver (flag: plot_network_with_numbering)
num_nodes = length(nodes);
num_edges = length(edges);

if (plot_network_with_numbering == true)
    figure(1), hold on, grid on
    %view(180, -90)
    for i = 1:num_nodes
        %plot3(nodes(i,1), nodes(i,2), nodes(i,3), 'k.', 'MarkerSize', 12)
        %text(nodes(i,1), nodes(i,2), nodes(i,3), num2str(i))
        plot(nodes(i,1), nodes(i,2), 'k.', 'MarkerSize', 12)
        text(nodes(i,1), nodes(i,2), num2str(i))
    end

    for j = 1:length(edges)
        %plot3([nodes(edges(j,1),1) nodes(edges(j,2),1)], [nodes(edges(j,1),2) nodes(edges(j,2),2)], [nodes(edges(j,1),3) nodes(edges(j,2),3)], 'k', 'LineWidth', 1)
        %text((nodes(edges(j,1),1) + nodes(edges(j,2),1))/2, (nodes(edges(j,1),2) + nodes(edges(j,2),2))/2, (nodes(edges(j,1),3) + nodes(edges(j,2),3))/2, num2str(j), 'Color', 'g', 'FontSize', 15)
        plot([nodes(edges(j,1),1) nodes(edges(j,2),1)], [nodes(edges(j,1),2) nodes(edges(j,2),2)], 'k', 'LineWidth', 1)
        text((nodes(edges(j,1),1) + nodes(edges(j,2),1))/2, (nodes(edges(j,1),2) + nodes(edges(j,2),2))/2, num2str(j), 'Color', 'g', 'FontSize', 15)
    end
end

%--------------------------------------------------------------------------
%% --- Flow solver
PBCs = [aorta_BC_nodes(1) PDA1; aorta_BC_nodes(2) PDA2];

for j = 1:NvISVs
    PBCs = [PBCs; vISV_BC_nodes(j) PPCV];
end

[num_BCs dum] = size(PBCs);

% Vessel conductance (um^4-s/kg)
G = zeros(num_edges, TS);

for t = 1:TS
    for v = 1:num_edges
        G(v,t) = (pi*D(v,t)^4)/(128*mu*edges(v,3));
    end
end

% Vessel flow (um^3/s)
% Multiply by 3.6 to get uL/hr
Q = zeros(num_edges, TS);

% Nodal pressure (kg/s^2/um, 1e-6 Pa)
P = zeros(num_nodes, TS);

% Nodal connectivy array
node_conn = edges(:, 1:2);

% Solve for pressure and flow for each time step
for t = 1:TS
    % Assemble the system of equations to solve for pressure
    A = zeros(num_nodes, num_nodes);
    b = zeros(num_nodes, 1);
    p = zeros(num_nodes,1);

    for n = 1:num_nodes
        for v = 1:num_edges

            if (node_conn(v,2) == n)
                A(n,n) = A(n,n) - G(v,t);
                A(n,node_conn(v,1)) = A(n,node_conn(v,1)) + G(v,t);
            end

            if (node_conn(v,1) == n)
                A(n,n) = A(n,n) - G(v,t);
                A(n,node_conn(v,2)) = A(n,node_conn(v,2)) + G(v,t);
            end
        end
    end

    % Apply the boundary conditions
    for bc = 1:num_BCs
        node = PBCs(bc,1);
        BC = PBCs(bc,2);

        avg_G = mean(G(:,t));

        A(node,:) = [zeros(1,node-1) -avg_G zeros(1, num_nodes-node)];
        b(node) = -avg_G*BC;
    end

    % Solve for nodal pressures (kg/s^2/um, 1e-6 Pa)
    p = A\b;
    
    P(:,t) = p;
end

% Calculate flow (um^3/s)
for t = 1:TS
    for v = 1:num_edges
        node1 = node_conn(v,1);
        node2 = node_conn(v,2);

        Q(v,t) = -G(v,t)*(P(node2,t) - P(node1,t));
    end
end

% Calculate wall shear stress (kg/s^2/um, 1e-6 Pa)
tau = zeros(num_edges, TS);

for t = 1:TS
    for v = 1:num_edges
        node1 = node_conn(v,1);
        node2 = node_conn(v,2);
        L = edges(v,3);
        
        tau(v,t) = (D(v,t)/(4*L))*abs(P(node2,t) - P(node1,t));
    end
end

% Convert pressure to Pa
P = P*1e6;

% Convert flow to uL/hr
Q = Q*3.6e-6;

% Convert flow to mL/hr
%Q = Q/1000;

% Convert wall shear stress to Pa
tau = tau*1e6;

%--------------------------------------------------------------------------
%% --- Plot the network
% Create the video file
vid_file_name = ['flow_in_idealised_ISV_network_vid'];
vidObj = VideoWriter(vid_file_name,'MPEG-4');
vidObj.FrameRate = 1;
vidObj.Quality = 100;
open(vidObj);

% Define info for the color scale
% flow_max = max(max(Q));
% flow_min = min(min(Q));
% flow_zero = 1e-5;

red_max = [1 0.5 0.5];
red_min = [0.9 0 0];
blue_max = [0.5 0.5 1];
blue_min = [0 0 0.9];
zero_gray = [0.8 0.8 0.8];

for t = 1:TS
    flow_max = max(Q(:,t));
    flow_min = min(Q(:,t));
    flow_zero = 1e-5;
    
    
    if(ishandle(2))
        close(2);
    end
    
    figure(2), hold on, 
    
    %box on, grid on

    vess_colors = zeros(num_edges,3);
    
    for v = 1:num_edges
        if (Q(v,t) > flow_zero)
            xi = Q(v,t)/flow_max;
            vess_colors(v,:) = (1 - xi)*red_min + xi*red_max;
        end

        if (Q(v,t) < -flow_zero)
            xi = Q(v,t)/flow_min;
            vess_colors(v,:) = (1 - xi)*blue_min + xi*blue_max;
        end

        if (abs(Q(v,t)) < flow_zero)
            vess_colors(v,:) = zero_gray;
        end
    end
        
    %plot3(nodes(:,1), nodes(:,2), nodes(:,3), 'k.', 'MarkerSize', 15)
    plot(nodes(:,1), nodes(:,2), 'k.', 'MarkerSize', 15)
    
    for j = 1:num_edges
%         x1 = nodes(edges(j,1),1);
%         y1 = nodes(edges(j,1),2);
%         z1 = nodes(edges(j,1),3);
%         x2 = nodes(edges(j,2),1);
%         y2 = nodes(edges(j,2),2);
%         z2 = nodes(edges(j,2),3);
%         
%         plot3([x1 x2], [y1 y2], [z1 z2], 'Color', vess_colors(j,:), 'LineWidth', D(j,t))
%         arrow3([x1 y1 z1], [x2 y2 z2], 'y', 0.8)
        
        x1 = nodes(edges(j,1),1);
        y1 = nodes(edges(j,1),2);
        x2 = nodes(edges(j,2),1);
        y2 = nodes(edges(j,2),2);
        
        plot([x1 x2], [y1 y2], 'Color', vess_colors(j,:), 'LineWidth', D(j,t))
        plot([x1 x2], [y1 y2], 'w.-', 'LineWidth', 1, 'MarkerSize', 10)
        
        %arrow3([x1 y1 0], [x2 y2 0], 'w', 0.6)
        
    end
    
    flow_range = [linspace(flow_min, -flow_zero, 4) linspace(flow_zero, flow_max, 4)]';
    
    leg_plot_length = (max(nodes(:,2)) - min(nodes(:,1)))/length(flow_range);
    leg_colormap = [];
    
    for i = 1:length(flow_range)
        if (flow_range(i) < -flow_zero)
            xi = flow_range(i)/flow_min;
            
            leg_plot_color = (1 - xi)*blue_min + xi*blue_max;
            
            %plot([x_max-x_range/10 x_max-x_range/10], [min(nodes(:,2))+(i-1)*leg_plot_length min(nodes(:,2))+i*leg_plot_length], 'Color', leg_plot_color, 'LineWidth', 10)
        end
        
        if (abs(flow_range(i)) <= flow_zero)
            leg_plot_color = zero_gray;
            
            %plot([x_max-x_range/10 x_max-x_range/10], [min(nodes(:,2))+(i-1)*leg_plot_length min(nodes(:,2))+i*leg_plot_length], 'Color', leg_plot_color, 'LineWidth', 10)
        end
        
        if (flow_range(i) > flow_zero)
            xi = flow_range(i)/flow_max;
            
            leg_plot_color = (1 - xi)*red_min + xi*red_max;
            
            %plot([x_max-x_range/10 x_max-x_range/10], [min(nodes(:,2))+(i-1)*leg_plot_length min(nodes(:,2))+i*leg_plot_length], 'Color', leg_plot_color, 'LineWidth', 10)
        end
        
        leg_colormap = [leg_colormap; leg_plot_color];
    end
    
    leg_tick_labels = cell(length(flow_range)+1,1);
    
    for i = 1:length(flow_range)/2
        leg_tick_labels{i} = num2str(flow_range(i));
    end
    
    leg_tick_labels{length(flow_range)/2 + 1} = num2str(0);
    
    for i = length(flow_range)/2 + 2:length(flow_range)+1
        leg_tick_labels{i} = num2str(flow_range(i-1));
    end
    
    colormap(leg_colormap);
    c = colorbar;  
    set(c, 'Ticks', linspace(0,1,length(flow_range)+1));
    set(c, 'TickLabels', leg_tick_labels);
    c.Label.String = ' Flow (\muL/hr)';
    c.Label.FontSize = 14;
    c.Label.FontName = 'Helvtica'
    
    %title([' Time Step ', num2str(t)])
    set(gcf, 'Color', 'w')
    set(gca, 'FontSize', 12)
    set(gca, 'LineWidth', 2)
    set(gca, 'XTickLabel', [])
    set(gca, 'YTickLabel', [])
    set(gca, 'ZTickLabel', [])
    
    axis off
    
    axis([min(nodes(:,1)) max(nodes(:,1)) min(nodes(:,2)) max(nodes(:,2))])
    
    fig2 = gcf;
    pos = fig2.Position;
    set(gcf, 'Position', [10 100 1200 400])
    
    % Write to the video file
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame)
    
    if (t == TS)
        i = 0;
        
        while (i < 3)
            % Write to the video file
            currFrame = getframe(gcf);
            writeVideo(vidObj,currFrame)
            i = i + 1;
        end
    end
    
    pause(0.2)
end

close(vidObj);

%% --- Create output plots
% Calculate aortic flow
Q_aorta_in = Q(1,:);

end_aorta_edge = 1;

for v = 1:num_edges
    if (vess_type(v) == "aorta")
        end_aorta_edge = v;
    end
end

Q_aorta_out = Q(end_aorta_edge,:);

Q_ISVs = Q_aorta_in - Q_aorta_out;
%Q_ISVs = Q_ISVs/(NaISVs + NvISVs);
   

% Pull out diameter, flow, and WSS for the ISVs
D_aISV = [];
Q_aISV = [];
tau_aISV = [];

for v = 1:num_edges
    if (vess_type(v) == "aISV")
        D_aISV = [D_aISV; D(v,:)];
        Q_aISV = [Q_aISV; Q(v,:)];
        tau_aISV = [tau_aISV; tau(v,:)];
    end
end

Q_aISV_end = [];
tau_aISV_end = [];

for v = 1:3:length(Q_aISV)
    Q_aISV_end = [Q_aISV_end; Q_aISV(v,:)];
    tau_aISV_end = [tau_aISV_end; tau_aISV(v,:)];
end

D_vISV = [];
Q_vISV = [];
tau_vISV = [];

for v = 1:num_edges
    if (vess_type(v) == "vISV")
        D_vISV = [D_vISV; D(v,:)];
        Q_vISV = [Q_vISV; Q(v,:)];
        tau_vISV = [tau_vISV; tau(v,:)]; 
    end
end

Q_vISV_end = [];
tau_vISV_end = [];

for v = 1:3:length(Q_vISV)
    Q_vISV_end = [Q_vISV_end; Q_vISV(v,:)];
    tau_vISV_end = [tau_vISV_end; tau_vISV(v,:)];
end

% --- Fig 3: Flow for all ISVs
figure(3), hold on, grid on, box on

aISV_num = linspace(1, NaISVs, NaISVs);
vISV_num = linspace(1, NvISVs, NvISVs);

plot(aISV_num, Q_aISV_end, 'r', 'LineWidth', 2)
plot(vISV_num, Q_vISV_end, 'b', 'LineWidth', 2)

axis([1 NaISVs -0.5 0.5])
xlabel(' ISV number ')
ylabel(' ISV flow (\muL/hr) ')
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 16)
set(gca, 'LineWidth', 2)

% -- Fig 4: Flow in middle 3 ISVs
figure(4), hold on, grid on, box on

plot(aISV_num(7:9), Q_aISV_end(7:9,:), 'r', 'LineWidth', 2)
plot(vISV_num(7:9), Q_vISV_end(7:9,:), 'b', 'LineWidth', 2)

title(' flow in last 3 ISVs ')
axis([7 9 -0.5 0.5])
xlabel(' ISV number ')
ylabel(' ISV flow (\muL/hr) ')
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 16)
set(gca, 'LineWidth', 2)

% -- Fig 5: Flow in middle 3 ISVs over time
figure(5), hold on, grid on, box on

plot(Q_aISV_end(7,:), 'r', 'LineWidth', 2)
plot(Q_aISV_end(8,:), 'r', 'LineWidth', 2)
plot(Q_aISV_end(9,:), 'r', 'LineWidth', 2)
plot(Q_vISV_end(7,:), 'b', 'LineWidth', 2)
plot(Q_vISV_end(8,:), 'b', 'LineWidth', 2)
plot(Q_vISV_end(9,:), 'b', 'LineWidth', 2)

title(' flow in middle 3 ISVs ')
axis([1 TS -0.1 0.1])
xlabel(' time steps ')
ylabel(' ISV flow (\muL/hr) ')
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 16)
set(gca, 'LineWidth', 2)

% --- Fig 6: WSS for all ISVs
figure(6), hold on, grid on, box on

aISV_num = linspace(1, NaISVs, NaISVs);
vISV_num = linspace(1, NvISVs, NvISVs);

plot(aISV_num, tau_aISV_end, 'r', 'LineWidth', 2)
plot(vISV_num, tau_vISV_end, 'b', 'LineWidth', 2)

axis([1 NaISVs 0 2])
xlabel(' ISV number ')
ylabel(' ISV WSS (Pa) ')
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 16)
set(gca, 'LineWidth', 2)

% -- Fig 7: WSS in middle 3 ISVs
figure(7), hold on, grid on, box on

plot(aISV_num(7:9), tau_aISV_end(7:9,:), 'r', 'LineWidth', 2)
plot(vISV_num(7:9), tau_vISV_end(7:9,:), 'b', 'LineWidth', 2)

title(' WSS in middle 3 ISVs ')
axis([7 9 0 2])
xlabel(' ISV number ')
ylabel(' ISV WSS (Pa) ')
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 16)
set(gca, 'LineWidth', 2)

% -- Fig 8: WSS in middle 3 ISVs over time
figure(8), hold on, grid on, box on

plot(tau_aISV_end(7,:), 'r', 'LineWidth', 2)
plot(tau_aISV_end(8,:), 'r', 'LineWidth', 2)
plot(tau_aISV_end(9,:), 'r', 'LineWidth', 2)
plot(tau_vISV_end(7,:), 'b', 'LineWidth', 2)
plot(tau_vISV_end(8,:), 'b', 'LineWidth', 2)
plot(tau_vISV_end(9,:), 'b', 'LineWidth', 2)

title(' WSS in middle 3 ISVs ')
axis([1 TS 0 2])
xlabel(' time steps ')
ylabel(' ISV WSS (Pa) ')
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 16)
set(gca, 'LineWidth', 2)

% -- Fig 9: ISV diameter
figure(9), hold on, grid on, box on

plot(DaISV, 'r', 'LineWidth', 2)
plot(DvISV, 'b', 'LineWidth', 2)

axis([1 TS 0 15])
xlabel(' time steps ')
ylabel(' ISV diameter (\mum) ')
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 16)
set(gca, 'LineWidth', 2)

Q_aISV_mid3 = [Q_aISV_end(7:9,:)]
Q_vISV_mid3 = [Q_vISV_end(7:9,:)]
tau_aISV_mid3 = [tau_aISV_end(7:9,:)]
tau_vISV_mid3 = [tau_vISV_end(7:9,:)]

