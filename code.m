% Soumya Ganguly
% PID : A53274333
% Final project 
% MAE 288 B: Optimal Estimation
clc;
clear all;

%% 
% Generation of data

% Number of timesteps
n_t = 1000;

% Noise parameters

prompt = 'Please input the range of the uniform distribution of disturbance';
d = input(prompt)

prompt = 'Please input the range of the uniform distribution of noise';
h=input(prompt)

Y_true = zeros(1, n_t);
X_true = zeros(1, n_t); 
X_true(1,1) = -pi + 2*pi*rand;
Y_true(1,1)=X_true(1, 1)^3+(-h + 2*h*rand);
for i = 2 : n_t
    X_true(1, i) = sin(X_true(1, i-1))+(-d + 2*d*rand); % generation of true state values
    Y_true(1, i) = X_true(1, i)^3+(-h + 2*h*rand); %generation of measurement
    
end
 
%%
% Bayesian filter

% Number of sample points for pdf
n_samp = 1000;

% Initializing filter
pdfX_f = zeros(n_t, n_samp);

X0 = linspace(-pi, pi, n_samp);
X = linspace(-1-d, 1+d, n_samp);

% Initializing predictor
pdfX_p = zeros(n_t, n_samp);
pdfX_p(1,:)= ones(1, n_samp)/(2*pi);

pdfX_f(1,:)= (1/(2*h))*(abs(Y_true(1, 1)-X0.^3) < h).*pdfX_p(1,:);
pdfX_f(1,:)= (n_samp/(2*pi))*pdfX_f(1,:)/sum(pdfX_f(1,:));
pdfX_p(2, :) = (2*pi/n_samp)*((1/(2*d))*(abs(X'-sin(X0)) < d)*pdfX_f(1,:)')';
    
X_mean_filtered(1,1)=((2*pi/n_samp)*pdfX_f(1,:)*X0');
var_f_X(1,1)=((2*pi/n_samp)*pdfX_f(1,:)*(X0.^2)') -(((2*pi/n_samp)*pdfX_f(1,:)*X0')^2);
X_mean_predicted(1,1)=((2*pi/n_samp)*pdfX_p(1,:)*X0');

for i = 2 : n_t
    % Filter update
    pdfX_f(i, :) = (1/(2*h))*(abs(Y_true(1, i)-X.^3) < h).*pdfX_p(i, :);
    pdfX_f(i, :) = (n_samp/(2*(1+d)))*pdfX_f(i, :)/(sum(pdfX_f(i, :)));
    % Predictor update
    pdfX_p(i+1, :) = (2*(1+d)/n_samp)*((1/(2*d))*(abs(X'-sin(X)) < d)*pdfX_f(i, :)')';
    X_mean_filtered_bayes(1,i) = ((2*(1+d)/n_samp)*pdfX_f(i,:)*X')';
    var_f_X(1,i)=((2*(1+d)/n_samp)*pdfX_f(i,:)*(X.^2)') -(((2*(1+d)/n_samp)*pdfX_f(i,:)*X')^2);
    X_mean_predicted_bayes(1,i)= ((2*(1+d)/n_samp)*pdfX_p(i,:)*X')';
end


figure(1)
plot(X_true)
hold on
plot(X_mean_predicted_bayes)
hold on
plot(X_mean_filtered_bayes)
legend('true values','predicted bayesian mean','filtered bayesian mean')
xlabel('discrete time points') 
ylabel('State means')
title('TASK 4.1:Comparison of predicted and filtered means in Bayesian Estimation')


figure(2)
plot((X_mean_predicted_bayes - X_true).^2)
hold on
plot((X_mean_filtered_bayes-X_true).^2)
legend('predicted bayesian mean squarred errors','filtered bayesian mean squarred error')
xlabel('discrete time points') 
ylabel('Mean Squarred Errors')
title('TASK 4.2:Comparison of predicted and filtered mean squarred errors')
%% Extended Kalman filter

% Initializing predictor
Pp = zeros(1, n_t+2);         %predictor covariance
xp = zeros(1, n_t+2);         % predictor state
yp = zeros(1, n_t+2);         % predictor output
Pp(1, 1) = pi^2/3;            % initial predictor covariance for uniform distribution
xp(1, 1) = 1e-1;              % initial predictor mean value, cannot take as zero, 
% because then every subsequent result will be zero

% Initializing filter
Pf = zeros(1, n_t+1); %filter covariance
xf = zeros(1, n_t+1); %filtered state
yf = zeros(1, n_t+1); % filtered output
% Variance of noise
Q = d^2/3; % true for uniform distribution
R = h^2/3;

% Run the Extended Kalman filter
for i = 1 : n_t
    % Filter update
    H = 3*xp(1, i)^2;
    K = Pp(1, i)*H/(H^2*Pp(1, i)+R); % kalman gain
    xf(1, i) = xp(1, i) + K*(Y_true(1, i)-xp(1, i)^3);
    Pf(1, i) = Pp(1, i)-Pp(1, i)^2*H^2/(H^2*Pp(1, i)+R);
    
    % Predictor update
    xp(1, i+1) = sin(xf(1, i));
    F = cos(xf(1, i));
    Pp(1, i+1) = F^2*Pf(1, i)+Q;
end

%Plots of Task 6

figure(3)
plot(X_true)
hold on
% plotting expected value of state obtained from bayesian filter
plot(X_mean_filtered_bayes)
hold on
%plotting EKF state data
plot(xf(1, 2 : n_t+1))
legend('true state values','bayesian state means','filtered values of Extended Kalman filter')
xlabel('discrete time points') 
ylabel('State mean/true values') 
title('TASK 6.1:Comparison of estimation of different methods with respect to the means and true state values')

figure(4)
plot(var_f_X)
%legend('variance of bayesian filter')
xlabel('discrete time points') 
ylabel('Variance') 
title('TASK 6.2:Variance of Bayesian Filter')

figure(5)
plot(Pf(1,2:n_t+1),'r')
xlabel('discrete time points') 
ylabel('Variance') 
%legend('variance of Extended Kalman Filter')
title('TASK 6.2:Variance of Extended Kalman Filter')
