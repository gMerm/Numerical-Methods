%example
n = 5;  
D = diag(4*ones(1, n));  
L = diag(-ones(1, n-1), -1);  
U = diag(-ones(1, n-1), 1);  
b = randn(n, 1);  

D_hat = 4 * ones(1, n);  
L_hat = -ones(1, n-1);  
U_hat = -ones(1, n-1);  

% Ensure that the sizes match
if length(L_hat) ~= n-1 || length(U_hat) ~= n-1
    error('Invalid sizes for L_hat or U_hat.');
end

%solve the tridiagonal system
x = solveTridiagonal(D, L, U, b);
y = solveTridiagonalHad(D_hat, L_hat, U_hat, b);

%res
disp('Solution x:');
disp(x);
disp('Solution y:');
disp(y);



%%
%larger problem
n = 1000;
D_large = diag(4 * ones(1, n));
L_large = diag(-ones(1, n-1), -1);
U_large = diag(-ones(1, n-1), 1);
b_large = randn(n, 1);

%non-vectorized function
tic;
x_non_vectorized = solveTridiagonal(D_large, L_large, U_large, b_large);
time_non_vectorized = toc;

%vectorized function
tic;
x_vectorized = solveTridiagonalHad(D_large, L_large, U_large, b_large);
time_vectorized = toc;

%execution times
disp(['Execution time (non-vectorized): ' num2str(time_non_vectorized) ' seconds']);
disp(['Execution time (vectorized): ' num2str(time_vectorized) ' seconds']);



%%
%diagonally dominant matrices C1-4
C1 = diag(10 * ones(1, n)) + diag(-3 * ones(1, n-1), 1) + diag(-3 * ones(1, n-1), -1);
C2 = diag(8 * (-1).^(1:n)) + diag(-2 * ones(1, n-1), 1) + diag(-2 * ones(1, n-1), -1);
C3 = diag(6 * ones(1, n)) + diag(-1.5 * (0.9).^(1:n-1), 1) + diag(-1.5 * (0.9).^(1:n-1), -1);
C4 = diag(7 * ones(1, n)) + diag(-2 * rand(1, n-1), 1) + diag(-2 * rand(1, n-1), -1);

x_vectorized_diag_dom1 = solveTridiagonalHad(diag(C1), diag(C1, -1), diag(C1, 1), b);
x_vectorized_diag_dom2 = solveTridiagonalHad(diag(C2), diag(C2, -1), diag(C2, 1), b);
x_vectorized_diag_dom3 = solveTridiagonalHad(diag(C3), diag(C3, -1), diag(C3, 1), b);
x_vectorized_diag_dom4 = solveTridiagonalHad(diag(C4), diag(C4, -1), diag(C4, 1), b);

disp('Solutions for C1,C2,C3,C4:');
disp([x_vectorized_diag_dom1,x_vectorized_diag_dom2,x_vectorized_diag_dom3,x_vectorized_diag_dom4]);

%calculate norms
norm_C1_upper = norm(triu(C1, 1)); 
norm_C1_lower = norm(tril(C1, -1));

norm_C2_upper = norm(triu(C2, 1)); 
norm_C2_lower = norm(tril(C2, -1)); 

norm_C3_upper = norm(triu(C3, 1)); 
norm_C3_lower = norm(tril(C3, -1)); 

norm_C4_upper = norm(triu(C4, 1)); 
norm_C4_lower = norm(tril(C4, -1)); 

%display norms
disp('Norms of Upper and Lower Triangular Parts:');
disp(['C1: Norm Upper = ', num2str(norm_C1_upper), ', Norm Lower = ', num2str(norm_C1_lower)]);
disp(['C2: Norm Upper = ', num2str(norm_C2_upper), ', Norm Lower = ', num2str(norm_C2_lower)]);
disp(['C3: Norm Upper = ', num2str(norm_C3_upper), ', Norm Lower = ', num2str(norm_C3_lower)]);
disp(['C4: Norm Upper = ', num2str(norm_C4_upper), ', Norm Lower = ', num2str(norm_C4_lower)]);

function x = solveTridiagonal(D, L, U, b)
    %D, L, U are the tridiagonal matrices
    %b is the right-hand side vector(s)

    n = length(b);  %size of the system
    x = zeros(n, size(b, 2));  %initialize the solution matrix

    %LU decomposition
    for k = 2:n
        factor = L(k, k-1) / D(k-1, k-1);
        D(k, k) = D(k, k) - factor * U(k-1, k);
        b(k, :) = b(k, :) - factor * b(k-1, :);
    end

    %backward substitution
    x(n, :) = b(n, :) / D(n, n);
    for k = n-1:-1:1
        x(k, :) = (b(k, :) - U(k, k+1) * x(k+1, :)) / D(k, k);
    end
end

function x = solveTridiagonalHad(D, L, U, b)
    % D, L, U are vectors of diagonals
    % b is the right-hand side matrix with one or more columns

    n = length(b);
    x = zeros(n, size(b, 2));  % Initialize the solution matrix

    %LU decomposition
    for k = 2:n
        factor = L(k-1) ./ D(k-1);  
        D(k) = D(k) - factor .* U(k-1);  
        b(k, :) = b(k, :) - factor .* b(k-1, :);
    end

    %Backward substitution
    x(n, :) = b(n, :) ./ D(n);  
    for k = n-1:-1:1
        x(k, :) = (b(k, :) - U(k) .* x(k+1, :)) ./ D(k);  
    end
end
