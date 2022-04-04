% CS754-2022 Assignment 4, Q1
% Arpon Basu and Shashwat Garg
%  200050013 and 200050130
addpath('./l1_ls_matlab');
rng(13)

n=500;
m=200;
p=0.5;
sqrt_m = sqrt(m);
x = zeros(n,1);
phi = 2*binornd(1,p,m,n)-1;
phi = phi/sqrt_m;

V = phi(1:m/10,:);
R = phi(1+m/10:m,:);

random_range = 1000;
for i=1:18
    x(randi(n))= unifrnd(0,random_range);
end
while nnz(x)<18 
    x(randi(n))= unifrnd(0,random_range);
end

% All thinfs initialised
y = phi*x;

noise_sigma = abs(0.05*sum(y,"all")/m);

y = y + normrnd(0,noise_sigma,m,1);
% Noise added to y

y_V = y(1:m/10,:);
y_R = y(1+m/10:m,:);

% Separated into R and V
lambda_set = [0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5 1 2 5 10 15 20 30 50 100];
estimated_x = zeros(length(lambda_set), n,1);
validation_error = zeros(length(lambda_set),1);
RMSE = zeros(length(lambda_set),1);

for i=1:length(lambda_set)
    [estimated_x(i,:,:),~]=l1_ls(R,y_R,lambda_set(i));
    RMSE(i) = norm(x-estimated_x(i,:,:)')/norm(x);
    validation_error(i) = sumsqr(y_V - V*estimated_x(i,:,:)')/length(y_V);
end
% Compute the results and errors for each value of lambda from the set.

format short
fprintf('\n\tLambda\t\tRMSE\tValidation Error\n')
disp([lambda_set' RMSE validation_error])
% Display the results

log_lambda = log(lambda_set);

figure()
fig1 = plot(log_lambda,RMSE);
title('RMSE wrt Log of Lambda')
xlabel('Log(Lambda)')
ylabel('RMSE')
saveas(fig1,'RMSE.png')

figure()
fig2=plot(log_lambda,validation_error);
title('VE wrt Log of Lambda')
xlabel('Log(Lambda)')
ylabel('Validation Error')
saveas(fig2,'VE.png')

% Plot the RMSE and VE wrt the logarithm of the lambda value.




