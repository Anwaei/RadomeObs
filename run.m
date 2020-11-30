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
close all
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

tol = 50;  % Min distance when guidance ends

% ------------ Filter settings ------------
dt = 0.05; t = 0:dt:25; T = dt;
num_steps = numel(t);
x_true = zeros(D,num_steps); x_true(:,1) = x0;  % Ture states
z = zeros(1,num_steps);  % Observations
ucmd = zeros(1,num_steps);  % Guidance input
u = zeros(1,num_steps);  % Control input
m = zeros(D,num_steps);  % Estimated means
P = zeros(D,D,num_steps);  % Estimated covariances
m(:,1) = [deg2rad(7.7106) 12050 deg2rad(2) 0]';  % Initial estimated mean
P(:,:,1) = diag([(deg2rad(1))^2 1000^2 (deg2rad(1))^2 1^2]);  % Initial estimated covariance
Q = zeros(D,D,num_steps);  % Process noise
R = zeros(1,num_steps);  % Observation noise
for k = 1:num_steps
    Q(:,:,k) = diag([(deg2rad(0.01))^2 10^2 (deg2rad(0.01))^2 0.001^2]);
    R(k) = (deg2rad(0.01))^2;
end
X_M = zeros(1,num_steps); X_M(1) = X_M0;
Y_M = zeros(1,num_steps); Y_M(1) = Y_M0;
V_Mx = zeros(1,num_steps); V_Mx(1) = V_Mx0;
V_My = zeros(1,num_steps); V_My(1) = V_My0;
lambda_dot = zeros(1,num_steps);
gamma_dot = zeros(1,num_steps);
u_int = zeros(1,num_steps);

% For test
error_R = zeros(1,num_steps);

% ------------ Model functions ------------
f = @dynFunc; h = @measFunc;
Fx = @dynJacob; Hx = @measJacob;

%% Simulation
for k = 1:num_steps
    % ----------- Filtering -----------
    if k > 1
        % ---------    EKF    ---------
        % Predict
        m_predict = f(m(:,k-1),u(:,k-1));
        P_predict = Fx(m(:,k-1),u(:,k-1))*P(:,:,k-1)*Fx(m(:,k-1),u(:,k-1))' + Q(:,:,k-1);
        % Updating
        noise_obs = mvnrnd(0,R(k));
        z_k = h(x_true(:,k)) + noise_obs;  % Current observation
        v_k = z_k - h(m_predict);  % Innovation
        S_k = Hx(m_predict) * P_predict * Hx(m_predict)' + R(k);
        K_k = P_predict * Hx(m_predict)' / S_k;  % Kalman gain
        m_updated = m_predict + K_k*v_k;
        P_updated = P_predict - K_k*S_k*K_k';
        
        z(k) = z_k;
        m(:,k) = m_updated;
        P(:,:,k) = P_updated;
    end
    
%     if abs(x_true(2,k)) < tol
    if norm([X_M(k)-X_T, Y_M(k)-Y_T]) < tol
        break;
    end
    V_Mx(k) = V_M*cos(x_true(3,k));
    V_My(k) = V_M*sin(x_true(3,k));
    
    % ---------- Propogating ----------
    lambda_dot(k) = -V_M/x_true(2,k)*sin(x_true(3,k)-x_true(1,k));
    ucmd(k) = N*V_M*lambda_dot(k);  % Proportional guidance
    if k < num_steps 
        u(:,k+1) = (tau-T)/tau * u(:,k) + T/tau*ucmd(:,k);  % Pilot input, using forward difference
        u_int(k+1) = u_int(k) + abs(u(:,k+1))*dt;
        noise_pro = mvnrnd(zeros(4,1),Q(:,:,k))';        
%         x_true(:,k+1) = f(x_true(:,k),u(:,k)) + noise_pro;  % State propagation
        x_true(:,k+1) = f(x_true(:,k),u(:,k));
        X_M(k+1) = X_M(k) + V_Mx(k)*dt;
        Y_M(k+1) = Y_M(k) + V_My(k)*dt;
    end
    gamma_dot(k) = u(:,k)/V_M;
    
    error_R(k) = norm([X_M(k)-X_T, Y_M(k)-Y_T]) - x_true(2,k);
    
end

%% Results
num_steps_end = k;
x_true = x_true(:,1:num_steps_end); 
X_M = X_M(:,1:num_steps_end); Y_M = Y_M(:,1:num_steps_end);
m = m(:,1:num_steps_end); P = P(:,:,1:num_steps_end);
t = t(:,1:num_steps_end); z = z(:,1:num_steps_end);
u = u(:,1:num_steps_end); ucmd = ucmd(:,1:num_steps_end); 
u_int = u_int(:,1:num_steps_end);
lambda_dot = lambda_dot(:,1:num_steps_end);
lambda_dot = rad2deg(lambda_dot);
gamma_dot = gamma_dot(:,1:num_steps_end);
gamma_dot = rad2deg(gamma_dot);
sigma_dot = gamma_dot - lambda_dot;
error_R = error_R(:,1:num_steps_end);
std = zeros(D, num_steps_end);
for k = 1:num_steps_end
    std(:,k) = sqrt(diag(P(:,:,k)));
end

zem = zeros(1,num_steps_end); t_end = t(end);  % Zero Effort Distance
for k = 1:num_steps_end
    t_go = t_end - t(k);
    X_end = X_M(k) + V_Mx(k)*t_go;
    Y_end = Y_M(k) + V_My(k)*t_go;
    zem(k) = sqrt((X_end-X_T)^2 + (Y_end-Y_T)^2);
end

x_true_deg = x_true;
x_true_deg(1,:) = rad2deg(x_true_deg(1,:));
x_true_deg(3,:) = rad2deg(x_true_deg(3,:));
m_deg = m;
m_deg(1,:) = rad2deg(m_deg(1,:));
m_deg(3,:) = rad2deg(m_deg(3,:));
std_deg = std;
std_deg(1,:) = rad2deg(std_deg(1,:));
std_deg(3,:) = rad2deg(std_deg(3,:));

state_labels = {'lambda', 'R', 'gamma_M', 'rho'};
for d = 1:D
    figure(d); hold on
    hi = patch([t, fliplr(t)],[m_deg(d,:)'-2*std_deg(d,:)'; flipud(m_deg(d,:)'+2*std_deg(d,:)')], 1, ...
        'FaceColor', [0.9,0.9,1], 'EdgeColor', 'none'); % Two times the standard deviation
    set(hi,'handlevisibility','off');
    plot(t,x_true_deg(d,:),t,m_deg(d,:));
    legend(['true ' state_labels{d}], ['estimated ' state_labels{d}])
    xlabel('t');
end
figure(D+1); hold on
plot(t,z); xlabel('t');
title('observations')
figure(D+2); hold on
plot(t,u,t,ucmd);
legend('input pilot', 'input guidance')
figure(D+3); hold on
plot(X_M,Y_M);
plot(X_T,Y_T,'o');
xlabel('x'); ylabel('y');
legend('missile','target')
title('trajectory')
figure(D+4); hold on
plot(t,lambda_dot);
title('lambda dot');

%% Figures in the paper
D1 = D+4;
figure(D1+1); hold on
plot(X_M,Y_M);
plot(X_T,Y_T,'o');
xlabel('x'); ylabel('y');
title('Missile Trajectory')
figure(D1+2); hold on
plot(t,lambda_dot);
title('Line-of-Sight Angle Rate');
figure(D1+3); hold on
plot(t,gamma_dot);
title('Missile Flight Path Angle Rate');
figure(D1+4); hold on
plot(t,sigma_dot);
title('Lead Angle Rate');
figure(D1+5); hold on
plot(t,zem);
title('Zero-Effort Miss Distance');
figure(D1+6); hold on
plot(t,u);
title('Missile Maneuver Acceleration');
figure(D1+7); hold on
plot(t,u_int);
title('Missile Maneuver Acceleration Integral');

