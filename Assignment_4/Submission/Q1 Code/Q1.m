% CS754-2022 Assignment 4, Q1
% Arpon Basu and Shashwat Garg

n=500;
m=200;
p=0.5;

x = zeros(n,1);
phi = 2*binornd(1,p,m,n)-1;

V = phi(1:m/10,:);
R = phi(1+m/10:m,:);

random_range = 1000;
for i=1:18
    x(randi(n))= unifrnd(0,random_range);
end
while nnz(x)<18 
    x(randi(n))= unifrnd(0,random_range);
end

y = phi*x;

noise_sigma = abs(0.05*sum(y,"all")/m);

y = y + normrnd(0,noise_sigma,m,1);

y_V = y(1:m/10,:);
y_R = y(1+m/10:m,:);
lambda_set = [0.0005 0.0001 0.005 0.001 0.01 0.05 0.1 0.5 1 2 5];
estimated_x = zeros(length(lambda_set), n,1);
validation_error = zeros(length(lambda_set),1);
RMSE = zeros(length(lambda_set),1);

for i=1:length(lambda_set)
    [estimated_x(i,:,:),~]=l1_ls(R,y_R,lambda_set(i));
    RMSE(i) = norm(x-estimated_x(i,:,:)')/norm(x);
    validation_error(i) = sumsqr(y_V - V*estimated_x(i,:,:)')/length(y_V);
end


format long
fprintf('\n\nRMSE\n')
disp(RMSE)
fprintf('Validation error\n')
disp(validation_error)









