tau = 0.1; dt = 0.01;
pilot_tf = tf(1,[tau,1]);
t1 = 0:0.01:2; % t_out = 0:0.005:1;
x1 = sin(5*t1);
y1 = lsim(pilot_tf,x1,t1);
t2 = 0:0.01:1;
x2 = sin(5*t2);
y2 = lsim(pilot_tf,x2,t2);
figure(1)
plot(t1,x1,t1,y1,t2,x2,t2,y2);
legend('x1','y1','x2','y2');

tau = 0.1; dt = 0.01;
t = 0:dt:2;
am = zeros(size(t)); ac = zeros(size(t));
for k = 1:numel(t)
    ac(k) = sin(5*t(k));
    if k < numel(t)
        am(k+1) = (tau-dt)/tau*am(k) + dt/tau*ac(k);
    end
%     if k > 1
%         am(k) = tau/(tau+dt)*am(k-1) + dt/(tau+dt)*ac(k);
%     end
end
figure(2)
plot(t,ac,t,am);
legend('ac','am')

K = 10000;
x = zeros(2,K);
L = [1,0.1;0.1,1];
for k = 1:K
    x(:,k) = mvnrnd(zeros(2,1),L);
end
plot(x(1,:),x(2,:),'.')
    