clear; clc;

%% ================== INITIALIZATION ==================
% Define material properties and geometric parameters with units
E = 200e9 * ones(7,1);  % Young's modulus (Pa or N/m^2) for all elements (Structural Steel)
A = 2 * ones(7,1);      % Cross-sectional area (m^2)
L = 8 * ones(7,1);      % Length of each truss member (m)
k = E.*A./L;            % Compute stiffness of each truss element (N/m)

% Define orientation angles (in radians) for each truss element
t = [0 0 pi/3 2*pi/3 pi/3 2*pi/3 0];

% Define connectivity of nodes (element-to-node mapping)
n = {[1 2], [2 3], [1 4], [2 4], [2 5], [3 5], [4 5]};

%% ============== VISUALIZATION ADDITION ==============
% Compute node positions for visualization
node_pos = zeros(5,2); % Initialize node positions (m)
node_pos(1,:) = [0 0]; % Set Node 1 at the origin

for elem = 1:length(n)
    nodes = n{elem};
    if all(node_pos(nodes(2),:) == 0) % If second node is not positioned
        theta = t(elem);
        length_m = L(elem); % Convert length to meters
        node_pos(nodes(2),:) = node_pos(nodes(1),:) + ...
                              length_m * [cos(theta), sin(theta)];
    end
end

% Plot the truss structure
figure('Name','Truss Visualization');
hold on; axis equal; grid on;
for elem = 1:length(n)
    nodes = n{elem};
    plot(node_pos(nodes,1), node_pos(nodes,2), 'b-o', 'LineWidth',1.5); % Plot truss elements
    text(mean(node_pos(nodes,1)), mean(node_pos(nodes,2)), num2str(elem), 'BackgroundColor','w'); % Label elements
end
scatter(node_pos(:,1), node_pos(:,2), 100, 'r', 'filled'); % Highlight nodes
text(node_pos(:,1)+0.5, node_pos(:,2), cellstr(num2str((1:5)'))); % Label nodes
hold off;

%% ============== TRUSS STIFFNESS COMPUTATION ==============
% Compute element stiffness matrices (N/m)
ke = truss_stiffness(k, t);  % Element stiffness matrices
K = truss_stiffness_g(k, t, n);  % Global stiffness matrix

%% ============== SOLVING FOR DISPLACEMENTS ==============
% Extract sub-matrix for reduced system of equations
% Partitioning the stiffness matrix to separate known and unknown displacements

rhs_forces = [0;0;0;0;0;0;-47880.25]; % Applied force vector (N)
K_reduced = [K(3:5,3:5) K(3:5,7:10); K(7:10,3:5) K(7:10,7:10)]; % Extract reduced stiffness matrix (N/m)
dl = K_reduced \ rhs_forces; % Solve for unknown displacements (m)

% Construct full displacement vector with known boundary conditions
d = [0;0;dl(1:3);0;dl(4:end)];

%% ============== COMPUTE FORCES AND STRESSES ==============
f = K * d;  % Compute nodal force vector (N)

% Compute stresses in each element (Pa or N/m^2)
s = zeros(length(E),1);
for cnt=1:length(s)
    f_n = 2 * n{cnt}(1) - 1; % X-index of first node
    s_n = 2 * n{cnt}(2) - 1; % X-index of second node
    
    % Compute stress using force balance in local coordinate system
    s(cnt,1) = (E(cnt) / L(cnt)) * [-1 1] * ...
              [cos(t(cnt)) sin(t(cnt)) 0 0;0 0 cos(t(cnt)) sin(t(cnt))] * ...
              [d(f_n:f_n+1);d(s_n:s_n+1)];
end

%% ============== DISPLAY RESULTS ==============

% Display global stiffness matrix (N/m)
fprintf('Total Stiffness Matrix (N/m):\n\n');
disp(K);

% Display nodal forces (N)
fprintf('\n--------\nNodal Forces (N):\n\n');
for cnt=1:length(f)/2
    fprintf('F%d_x = %7.3f N\n', cnt, f(2*cnt-1));
    fprintf('F%d_y = %7.3f N\n', cnt, f(2*cnt));
end

% Display nodal displacements (m)
fprintf('\n--------\nNodal Displacements (m):\n\n');
for cnt=1:length(d)/2
    fprintf('d%d_x = %7.3e m\n', cnt, d(2*cnt-1));
    fprintf('d%d_y = %7.3e m\n', cnt, d(2*cnt));
end

% Display element stresses (Pa or N/m^2)
fprintf('\n====== Element Stresses (Pa) ======\n');
for cnt=1:length(s)
    fprintf('Element %d: %+7.3e Pa\n', cnt, s(cnt));
end
