%question 1
%Measure execution times
n_values = 100:100:2000;
execution_times = zeros(size(n_values));

for i = 1:length(n_values)
    n = n_values(i);
    A = randn(n, n);
    A = A * A';  
    A = A + n * eye(n);
    timing = timeit(@() chol(A));
    execution_times(i) = timing;
end

%Fit a cubic function using polyfit
degree = 3; 
coefficients = polyfit(n_values, execution_times, degree);

%coefficients
disp('Cubic function coefficients:');
disp(coefficients);

%plot data and fitted curve   
fit_curve = polyval(coefficients, n_values);
figure;
plot(n_values, execution_times, 'o', n_values, fit_curve, '-');
xlabel('Matrix Size (n)');
ylabel('Execution Time');
title('Fitted Cubic Function for Execution Times using cubic polynomials');
legend('Measured Times', 'Fitted Cubic Curve');


%%
%question 2
%coefficients obtained from polyfit
coefficients = [0.0000; -0.0000; 0.0072; -0.9224] * 1e-03;

%values of n for prediction
n_values = 100:100:2000;
n_values1= 150:100:1550;

%calculate predicted execution times
predicted_times1 = polyval(coefficients, n_values);
predicted_times2 = polyval(coefficients, n_values1);

disp('Predicted Execution Times for 100:100:2000 (using cubic polynomials):');
disp(predicted_times1);
disp("");
disp('Predicted Execution Times for 150:100:1550 (using cubic polynomials):');
disp(predicted_times2);


%%
%question 3
%predicted execution times for both sets of n values
predicted_n_values = [100:100:2000, 150:100:1550];
predicted_execution_times = polyval(coefficients, predicted_n_values);

%visualize the results
figure;
plot(n_values, execution_times, 'o', predicted_n_values, predicted_execution_times, '-');
xlabel('Matrix Size (n)');
ylabel('Execution Time');
title('Actual vs Predicted Execution Times (using cubic polynomials)');
legend('Actual Times', 'Predicted Times', 'Location', 'NorthWest');



%%
%question 4, degree = 2
%Fit a 2d degree function using polyfit
degree1 = 2; 
coefficients1 = polyfit(n_values, execution_times, degree1);

%coefficients
disp('Function coefficients for 2nd degree polynomials:');
disp(coefficients1);

%plot data and fitted curve   
fit_curve1 = polyval(coefficients1, n_values);
figure;
plot(n_values, execution_times, 'o', n_values, fit_curve1, '-');
xlabel('Matrix Size (n)');
ylabel('Execution Time');
title('Fitted 2d Degree Polynomial for Execution Times');
legend('Measured Times', 'Fitted 2nd degree Curve');



%%
%question 4, degree = 4
%Fit a 4th degree function using polyfit
degree2 = 4; 
coefficients2 = polyfit(n_values, execution_times, degree2);

%coefficients
disp('Function coefficients for 4th degree polynomials:');
disp(coefficients2);

%plot data and fitted curve   
fit_curve2 = polyval(coefficients2, n_values);
figure;
plot(n_values, execution_times, 'o', n_values, fit_curve2, '-');
xlabel('Matrix Size (n)');
ylabel('Execution Time');
title('Fitted 4th Degree Polynomial for Execution Times');
legend('Measured Times', 'Fitted 4th degree Curve');
