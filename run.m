%% Model Configuration
% ----------------------------------------
% states: x = [lambda, R, gamma_M rho_theta]'
% input: u = A_M
% measurement: lambda_ME = (1+rho_theta)*lambda - rho_theta*gamma_M
% other states: X_M Y_M V_Mx V_My X_T Y_T
% guidance law: PNG: A_cmd = N*lambda_dot = -N/R*V_M*sin(gamma_M-lambda)
% discreted dynamic function: x(k+1) = f(x(k),u(k)) Jacobian function: Fx(x,u)
% discreted measurement function: z(k) = h(x(k)) Jacobian function: Hx(x)
% discreted pilot: A_M(k+1) = g(A_M(k),A_cmd(k))
% Stationary target
% Constant velocity
% Open loop: The guidance command is generated from the actual state variables
% ----------------------------------------

clear
global V_M T tau
% ------------ Initial states -------------
X_M0 = 0; Y_M0 = 0; V_Mx0 = 500; V_My0 = 0;
gamma_M0 = 0; X_T = 10^4; Y_T = 10^3;
rho_theta0 = 0.025;  % radome slope, assumed linear
V_M = norm([V_Mx0, V_My0]);

lambda0 = atan(Y_T/X_T);
R0 = norm([X_T, Y_T]);

tau = 0.1; N = 4;
x0 = [lambda0, R0, gamma_M0, rho_theta0]'; D = numel(x0);

% ------------ Filter settings ------------
dt = 0.05; t = 0:dt:25; T = dt;
num_steps = numel(t);
x_true = zeros(D,num_steps); x_true(:,1) = x0;
m = zeros(D,num_steps);
P = zeros(D,D,num_steps);
m(:,1) = [deg2rad(7.7106) 12050 deg2rad(2) 0]';  % Initial estimated mean
P(:,:,1) = diag([(deg2rad(1))^2 1000^2 (deg2rad(1))^2 1^2]);  % Initial estimated covariance
Q = zeros(D,D,num_steps);
R = zeros(1,num_steps);
for k = 1:num_steps
    Q(:,:,k) = diag([(deg2rad(0.01))^2 10^2 (deg2rad(0.01))^2 0.001^2]);
    R(k) = (deg2rad(0.01))^2;
end

% ------------ Model functions ------------
f = @dynFunc; h = @measFunc;
Fx = @dynJacob; Hx = @measJacob;

%% Simulation
