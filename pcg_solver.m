%load matrix from the path file
rng(42);
matrix_path = '1138_bus.mat';  
loaded_data = load(matrix_path);
A = loaded_data.Problem.A;

%random b vector for the right hand side
b = randn(size(A, 1), 1);

%tolerance & max iterations & sol
tol = 1e-6;
maxit =4;
xsol_ref=A\b;
M1 = ichol(A);
M3=ilu(A);
tic;
[x, flag, relres, iter, resvec, errvec] = pcg_1084639(A, b, tol, maxit,  [],[],[], xsol_ref);
solve_time = toc;


%display results
%disp('Solution x:');
%disp(x);
disp(['Flag: ', num2str(flag)]);
disp(['Relative Residual: ', num2str(relres)]);
disp(['Number of Iterations: ', num2str(iter)]);
disp(['Solving Time: ', num2str(solve_time), ' seconds']);


%calculate relative residual and error
rel_residual = resvec / norm(b);
rel_error = errvec / norm(xsol_ref);

%{
disp("Residual");
disp(rel_residual);
disp("Error");
disp(rel_error);
%}

%plot the convergence
figure;

%plot the relative residual
subplot(2, 1, 1);
semilogy(1:length(rel_residual), rel_residual, '-o');
title('Convergence Plot');
xlabel('Iteration');
ylabel('Relative Residual Norm (log scale)');
grid on;

%plot the relative error
subplot(2, 1, 2);
semilogy(1:length(rel_error), rel_error, '-o');
xlabel('Iteration');
ylabel('Relative Error Norm (log scale)');
grid on;

sgtitle('PCG Convergence Analysis');
disp(newline);

%%
%POISSON
n = 100;  %Number of internal nodes in each direction
m = 100;  %Number of points in the other direction


%Discretization
h = 1 / (n + 1);
x = linspace(0, 1, n + 2);
y = linspace(0, 0.5, m + 2);

%Create the matrix A
AA = sparse(n * m, n * m);

%loop over internal nodes (n)
for i = 2:n + 1
    %loop over points (m)
    for j = 2:m + 1
        %indices
        k = (j - 2) * n + i - 1;  

        %diagonal entry
        AA(k, k) = -4 / h^2;

        %neighbors
        if i > 2
            AA(k, k - 1) = 1 / h^2;
        end
        if i < n + 1
            AA(k, k + 1) = 1 / h^2;
        end
        if j > 2
            AA(k, k - n) = 1 / h^2;
        end
        if j < m + 1
            AA(k, k + n) = 1 / h^2;
        end
    end
end


bb = randn(size(AA, 1), 1);
xsol1 = AA \ bb;
M2=ilu(AA);

tic;
[x1, flag1, relres1, iter1, resvec1, errvec1] = pcg_1084639(AA, bb, tol, maxit, M2,[],[], xsol1);
solve_time1 = toc;


%display results
%disp('Solution x:');
%disp(x1);


disp(['Flag: ', num2str(flag1)]);
disp(['Relative Residual: ', num2str(relres1)]);
disp(['Number of Iterations: ', num2str(iter1)]);
disp(['Solving Time: ', num2str(solve_time1), ' seconds']);

%calculate relative residual and error
magnitude_errvec1 = abs(errvec1);
rel_residual_poisson = resvec1 / norm(bb);
rel_error_poisson = magnitude_errvec1 / norm(xsol1);

%{
disp("Residual");
disp(rel_residual_poisson);
disp("Error");
disp(rel_error_poisson);
%}

%plot the convergence
figure;

%plot the residual
subplot(2, 1, 1);
semilogy(1:length(rel_residual_poisson), rel_residual_poisson, '-o');
title('Convergence Plot');
xlabel('Iteration');
ylabel('Relative Residual Norm (log scale)');
grid on;

%plot the error
subplot(2, 1, 2);
semilogy(1:length(rel_error_poisson), rel_error_poisson, '-o');
xlabel('Iteration');
ylabel('Relative Error Norm (log scale)');
grid on;
sgtitle('PCG Convergence Analysis');


%%
function [x,flag,relres,iter,resvec,errvec] = pcg_1084639(A,b,tol,maxit,M1,M2,x0,xsol,varargin)
    
    if nargin < 2
        error(message('MATLAB:pcg:NotEnoughInputs'));
    end
    
    %MERM: Determine whether A is a matrix or a function.
    if isnumeric(A)
        atype = 'matrix';
        afun = @(x) A * x; 
    else
        atype = 'function_handle';
        afun = A;
    end

    if strcmp(atype,'matrix')
        %check matrix and right hand side vector inputs have appropriate sizes
        [m,n] = size(A);
        if (m ~= n)
            error(message('MATLAB:pcg:NonSquareMatrix'));
        end
        if ~isequal(size(b),[m,1])
            error(message('MATLAB:pcg:RSHsizeMatchCoeffMatrix', m));
        end
    else
        m = size(b,1);
        n = m;
        if ~iscolumn(b)
            error(message('MATLAB:pcg:RSHnotColumn'));
        end
    end
    
    %assign default values to unspecified parameters
    if (nargin < 3) || isempty(tol)
        tol = 1e-6;
    end
    warned = 0;
    if tol <= eps
        warning(message('MATLAB:pcg:tooSmallTolerance'));
        warned = 1;
        tol = eps;
    elseif tol >= 1
        warning(message('MATLAB:pcg:tooBigTolerance'));
        warned = 1;
        tol = 1-eps;
    end
    if (nargin < 4) || isempty(maxit)
        maxit = min(n,20);
    end
    maxit = max(maxit, 0);
    
    %check for all zero right hand side vector => all zero solution
    n2b = norm(b);                     %Norm of rhs vector, b
    if (n2b == 0)                      %if    rhs vector is all zeros
        x = zeros(n,1);                %then  solution is all zeros
        flag = 0;                      %a valid solution has been obtained
        relres = 0;                    %the relative residual is actually 0/0
        iter = 0;                      %no iterations need be performed
        resvec = 0;                    %resvec(1) = norm(b-A*x) = norm(0)
        if (nargout < 2)
            itermsg('pcg',tol,maxit,0,flag,iter,NaN);
        end
        return
    end
    
    %MERM
    if (nargin >= 5) && ~isempty(M1)
        existM1 = 1;
        if isnumeric(M1)
            m1type = 'matrix';
            m1fun = M1;
        elseif isa(M1, 'function_handle')
            m1type = 'function_handle';
            m1fun = M1;
        else
            error('Invalid preconditioner M1');
        end
    else
        existM1 = 0;
        m1type = 'matrix';
        m1fun = [];
    end
    
    %MERM
    if ((nargin >= 6) && ~isempty(M2))
        existM2 = 1;
        if isnumeric(M2)
            m2type = 'matrix';
            m2fun = M2;
        elseif isa(M2, 'function_handle')
            m2type = 'function_handle';
            m2fun = M2;
        else
            error('Invalid preconditioner M2');
        end
    else
        existM2 = 0;
        m2type = 'matrix';
    end
    
    if ((nargin >= 7) && ~isempty(x0))
        if ~isequal(size(x0),[n,1])
            error(message('MATLAB:pcg:WrongInitGuessSize', n));
        else
            x = x0;
        end
    else
        x = zeros(n,1);
    end
    
    %MERM
    if ((nargin > 9) && strcmp(atype,'matrix') && ...
            strcmp(m1type,'matrix') && strcmp(m2type,'matrix'))
        error(message('MATLAB:pcg:TooManyInputs'));
    end

    %MERM 
    if ((nargin >= 8) && ~isempty(xsol))
        if ~isequal(size(xsol), [n, 1])
            error(message('MATLAB:pcg:WrongExactSolSize', n));
        end
    else
        %empty xsol if not provided by the user
        xsol = []; 
    end
    
    %set up for the method
    flag = 1;
    xmin = x;                          %Iterate which has minimal residual so far
    imin = 0;                          %Iteration at which xmin was computed
    tolb = tol * n2b;                  %Relative tolerance
    
    %MERM
    %r = b - iterapp('mtimes',afun,atype,x,varargin{:});
    if strcmp(atype, 'matrix')
        
        r = b - A * x; 
    elseif strcmp(atype, 'function_handle')
        r = b - afun(x); 
    else
        error('Unsupported atype for matrix-vector multiplication.');
    end

    normr = norm(r);                   %Norm of residual
    normr_act = normr;
    
    if (normr <= tolb)                 %Initial guess is a good enough solution
        flag = 0;
        relres = normr / n2b;
        iter = 0;
        resvec = normr;
        if (nargout < 2)
            itermsg('pcg',tol,maxit,0,flag,iter,relres);
        end
        return
    end
    
    resvec = zeros(maxit+1,1);         %Preallocate vector for norm of residuals
    resvec(1,:) = normr;               %resvec(1) = norm(b-A*x0)
    normrmin = normr;                  %Norm of minimum residual
    rho = 1;
    stag = 0;                          %stagnation of the method
    moresteps = 0;
    maxmsteps = min([floor(n/50),5,n-maxit]);
    maxstagsteps = 3;
    
    

    %MERM: Initialize error vector
    errvec = zeros(maxit+1, 1); 
    
    for ii = 1 : maxit
        if existM1
            
            y = m1fun \ r;
            if ~allfinite(y)
                flag = 2;
                break
            end
        else 
            y = r;
        end
        
        if existM2
            
            z = m2fun \r;
            if ~allfinite(z)
                flag = 2;
                break
            end
        else 
            z = y;
        end
        
        rho1 = rho;
        rho = r' * z;
        if ((rho == 0) || isinf(rho))
            
            flag = 4;
            break
        end
        if (ii == 1)
            p = z;
        else
            beta = rho / rho1;
            if ((beta == 0) || isinf(beta))
                
                flag = 4;
                break
            end
            p = z + beta * p;
        end
        if strcmp(atype, 'matrix')
            
            q = A * p; 
        elseif strcmp(atype, 'function_handle')
            q = afun(p); 
        else
            error('Unsupported atype for matrix-vector multiplication.');
        end
        pq = p' * q;
        if (isinf(pq))
            
            flag = 4;
            break
        else
            alpha = rho / pq;
        end
        if isinf(alpha)
            
            flag = 4;
            break
        end
        
        %check for stagnation of the method    
        if (norm(p)*abs(alpha) < eps*norm(x))
            stag = stag + 1;
        else
            stag = 0;
        end
        
        x = x + alpha * p;             %form new iterate
        r = r - alpha * q;
        normr = norm(r);
        normr_act = normr;
        resvec(ii+1,1) = normr;
        
        %MERM: Calculate A-norm of the error from xsol to x
        if ~isempty(xsol)
            %Calculate the norm of the difference
            %errvec(ii+1, 1) = norm(x - xsol); 
            errvec(ii+1, 1) = sqrt((x - xsol)' * A * (x - xsol));
            
        else
            %If xsol is not provided, set error to NaN
            errvec(ii+1, 1) = NaN; 
        end

        %check for convergence
        if (normr <= tolb || stag >= maxstagsteps || moresteps)
            
            %MERM
            %r = b - iterapp('mtimes',afun,atype,x,varargin{:});
            if strcmp(atype, 'matrix')
                r = b - A * x;
            elseif strcmp(atype, 'function_handle')
                r = b - afun(x); 
            else
                error('Unsupported atype for matrix-vector multiplication.');
            end


            normr_act = norm(r);
            resvec(ii+1,1) = normr_act;

            if (normr_act <= tolb)
                flag = 0;
                iter = ii;
                break
            else
                if stag >= maxstagsteps && moresteps == 0
                    stag = 0;
                end
                moresteps = moresteps + 1;
                if moresteps >= maxmsteps
                    if ~warned
                        warning(message('MATLAB:pcg:tooSmallTolerance'));
                    end
                    flag = 3;
                    iter = ii;
                    break;
                end
            end

        end
        if (normr_act < normrmin)      %update minimal norm quantities
            normrmin = normr_act;
            xmin = x;
            imin = ii;
        end
        if stag >= maxstagsteps
            flag = 3;
            break;
        end
    end                                %for ii = 1 : maxit
    
    if isempty(ii)
        ii = 0;
    end
    
    %returned solution is first with minimal residual
    if (flag == 0)
        relres = normr_act / n2b;
    else

        %MERM
        %r_comp = b - iterapp('mtimes',afun,atype,xmin,varargin{:});
        if strcmp(atype, 'matrix')
            r_comp = b - A * xmin; 
        elseif strcmp(atype, 'function_handle')
            r_comp = b - afun(xmin); 
        else
            error('Unsupported atype for matrix-vector multiplication.');
        end

        if norm(r_comp) <= normr_act
            x = xmin;
            iter = imin;
            relres = norm(r_comp) / n2b;
        else
            iter = ii;
            relres = normr_act / n2b;
        end
    end
    

    %MERM: truncate the zeros from resvec and errvec
    %elements from the first row up to the current iteration: ii+1
    if ((flag <= 1) || (flag == 3))
        
        resvec = resvec(1:ii+1, :);
        errvec = errvec(1:ii+1, :);
        
    else
        
        resvec = resvec(1:ii, :);
        errvec = errvec(1:ii, :);
    end
    
    %only display a message if the output flag is not used
    if (nargout < 2)
        itermsg('pcg',tol,maxit,ii,flag,iter,relres);
    end

end
