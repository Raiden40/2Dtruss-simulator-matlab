function result = truss_stiffness_g(k,t,n)
% Function to compute the global stiffness matrix for a truss structure
% 
% Inputs:
% k - Vector of element stiffness values (N/m)
% t - Vector of element angles (radians)
% n - Cell array containing node connectivity for each element
% 
% Output:
% result - Global stiffness matrix (N/m)

% Ensure the stiffness and angle arrays have the same length
if length(k) ~= length(t)
    error('The two matrices must be of the same length.');
end

num_elem = length(k); % Determine the number of elements

% ---- Determine the total number of nodes ----
node_num = max(cellfun(@max, n)); % Find the highest node number in connectivity list

% Initialize the global stiffness matrix with zeros
K = zeros(2*node_num, 2*node_num);

% Loop through each truss element to assemble the global stiffness matrix
for cnt=1:num_elem
    C = cos(t(cnt)); % Compute cosine of the element angle
    S = sin(t(cnt)); % Compute sine of the element angle
    
    % Compute the element stiffness matrix in global coordinates
    kb = k(cnt) * [C^2   C*S   -C^2  -C*S;
                   C*S   S^2   -C*S  -S^2;
                  -C^2  -C*S    C^2   C*S;
                  -C*S  -S^2    C*S   S^2];

    % Identify the nodes for the current element
    if iscell(n)
        nodes = n{cnt}; % If n is a cell array
    else
        nodes = n(cnt,:); % If n is a numeric array
    end

    % Determine the global DOF indices for the two nodes of the element
    f_n = 2 * nodes(1) - 1; % First node index
    s_n = 2 * nodes(2) - 1; % Second node index
    
    % Assemble the element stiffness matrix into the global stiffness matrix
    K(f_n:f_n+1, f_n:f_n+1) = K(f_n:f_n+1, f_n:f_n+1) + kb(1:2,1:2);
    K(f_n:f_n+1, s_n:s_n+1) = K(f_n:f_n+1, s_n:s_n+1) + kb(1:2,3:4);
    K(s_n:s_n+1, f_n:f_n+1) = K(s_n:s_n+1, f_n:f_n+1) + kb(3:4,1:2);
    K(s_n:s_n+1, s_n:s_n+1) = K(s_n:s_n+1, s_n:s_n+1) + kb(3:4,3:4);
end

% Return the computed global stiffness matrix
result = K;

end