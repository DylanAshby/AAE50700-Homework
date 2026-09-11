% Define the matricies for the problem
A = [0 1 0; -1 0 0; 0 0 1];
B = [1 0 0; 0 1 0; 0 0 -1];
C = [0 1 0; -1 0.1 0; 0 0 1];

% Verify they are orthogonal
A*A'
B*B'
C*C'

% Verify they are not singular
det(A)
det(B)
det(C)