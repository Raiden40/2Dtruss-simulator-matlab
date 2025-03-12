function result = truss_stiffness(k, t)
    % This function computes the global stiffness matrix for each truss element.
    % Inputs:
    %   k - Vector of stiffness values for each element.
    %   t - Vector of angles (in radians) corresponding to each element.
    % Output:
    %   result - A 3D array containing the global stiffness matrices for all elements.

    % Check if the input vectors k and t have the same length
    if length(k) ~= length(t)
        error('The two matrices must be of the same length.');
    end

    num_elem = length(k); % Number of elements in the truss

    % Initialize a 3D matrix to store the global stiffness matrices
    % Each element's stiffness matrix is 4x4, and we store it in a 3D array.
    K = zeros(4, 4, num_elem);

    % Loop through each element to compute its global stiffness matrix
    for cnt = 1:num_elem
        C = cos(t(cnt)); % Cosine of the angle for current element
        S = sin(t(cnt)); % Sine of the angle for current element

        % Compute the global stiffness matrix for the current element
        kb = k(cnt) * [C^2  C*S  -C^2  -C*S;
                       C*S  S^2  -C*S  -S^2;
                      -C^2 -C*S   C^2   C*S;
                      -C*S -S^2   C*S   S^2];

        % Store the stiffness matrix in the 3D array
        K(:, :, cnt) = kb;
    end

    result = K; % Return the computed stiffness matrices
end
