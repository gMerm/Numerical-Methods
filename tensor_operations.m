%create and initialize an MDA tensor & a vector
rng(42);
tensor_MDA = randi([1, 10], 3, 3, 3);
vector = [2; 3; 4];
mode = 3;
%MDA tensors for ex2
A = [1, 2; 3, 4];
B = [5, 6; 7, 8];


test_tensor(tensor_MDA, vector, mode, A, B);

disp(newline);
disp("Given Test:")
test_tensor_given();





function result = ttv_1084639(tensor, vec, mode)
    disp(tensor);
    disp(vec);
    [m, n, p] = size(tensor);
    
    %check if matrix & vector are of the correct type
    if ~ismatrix(tensor) && ndims(tensor) > 2
        disp('Input tensor is a higher-dimensional array.');
    else
        error('Tensor must be a matrix or a higher-dimensional array.');
    end
    if ~isvector(vec)
        error('Vector must be a vector.');
    end

    %check if mode is 1,2 or 3
    if mode ~= 1 && mode ~= 2 && mode ~= 3
        error('Invalid mode. Mode must be 1, 2, or 3.');
    end
    
    %check if dims are compatible
    if numel(vec) ~= n
        error('Vector dimension does not match the size of the specified mode in the tensor.');
    end

    
    %initialize the result matrix
    result = zeros(m, p);
    
    if mode == 2
        %ttv operation for mode 2
        for i = 1:p
            result(:, i) = tensor(:,:,i) * vec;
        end
    elseif mode == 1
        %ttv operation for mode 1
        for i = 1:p
            result(:, i) = tensor(:,:,i).' * vec;
        end
    elseif mode == 3
        %ttv operation for mode 3
        for i = 1:m
            result(i, :) = tensor(i, :, 1) * vec(1) + tensor(i, :, 2) * vec(2) + tensor(i, :, 3) * vec(3);
        end
    end
end




function C = ttt_1084639(A, B, operation)
    %check the number of input arguments
    if nargin < 3
        operation = 'outer';
    end
    if nargin < 2
        error('Arguments number can''t be < 2.');
    end


    %check if input MDA tensors have the same number of dimensions
    if ndims(A) ~= ndims(B)
        error('Input tensors must have the same number of dimensions.');
    end
    
    %outer product that will be called if operations="outer" or if the
    %nargin < 3
    if strcmp(operation, 'outer')
        
        C = zeros([size(A), size(B)]);
        
        %in each iteration, the corresponding element of C is computed by
        %multiplying the element of index (i,j) from A and index (k,l)
        %from B. The result is a tensor with dims [size(A), size(B)].
        % i: rows of A, j: cols for A, k: rows of B, l: cols for B
        for i = 1:size(A, 1)
            for j = 1:size(A, 2)
                for k = 1:size(B, 1)
                    for l = 1:size(B, 2)
                        C(i, j, k, l) = A(i, j) * B(k, l);
                    end
                end
            end
        end
    
    %inner product
    elseif strcmp(operation, 'all')
        
        C = sum(sum(A .* B));

    else
        error('Invalid operation. Use ''outer'' or ''all''.');
    end
end


function test_tensor(tensorMDA, vector, mode, A, B)
    disp("TTV OPERATIONS:");
    %ttv_1084639 example
    disp("Result using my function and mode: "+ num2str(mode));

    %result using my function
    result = ttv_1084639(tensorMDA, vector, mode);
    disp(result);
   
    %"tensorize" the MDA tensor, so it can be used in TTs function ttv
    tensor_in=tensor(tensorMDA);
    tt_result = ttv(tensor_in, {vector},mode);
    disp("Result using Tensor's Toolbox function and mode: "+ num2str(mode));
    disp(tt_result);



    %ttt_1084639 example
    disp(newline);
    disp(newline);
    disp("PRODUCTS:");

    %outer prod
    result_outer = ttt_1084639(A, B);
    disp('Outer product result using ttt_1084639:')
    disp(result_outer)
    
    %inner prod
    result_inner = ttt_1084639(A, B, 'all');
    disp('Inner product result  using ttt_1084639:')
    disp(result_inner)
    
    %ttt using TensorToolbox
    tensorA=tensor(A);
    tensorB=tensor(B);
    result_ttt_outer=ttt(tensorA,tensorB);
    disp('Outer product result using Tensor''s Toolbox ttt function:')
    disp(result_ttt_outer);
    result_ttt_inner=ttt(tensorA,tensorB, [1,2]);
    disp('Inner product result using Tensor''s Toolbox ttt function:')
    disp(result_ttt_inner);
    
end


function test_tensor_given()
    clear; tol = 1e-8; nd = 3; rng(0); err = zeros(1, nd + 2);
    ndim = [2, 3, 4]; Atemp = randi(5, ndim); Btemp = randi(4, ndim); X = randi([-1, 1], max(ndim), 1); A = tensor(Atemp); B = tensor(Btemp);
    try
        for k = 1:nd
            err(k) = norm(ttv_1084639(A, X(1:ndim(k), 1), k) - double(ttv(A, X(1:ndim(k), 1), k)));
        end
        assert(max(err) < tol, 'ttv modal multiplication fails')
    catch ME1
    end
end
 


