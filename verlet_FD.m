function [u,e,error_mat,t] = verlet_FD(I,w,dt,T)
%%Solves the equation u_tt = -w^2 u, u(0) = I, u'(0) = 0 the standard vibration ODE using the
%%stromer method or verlet method. I and w are constants of any value, dt
%%is step size and T is the final time
t = [0:dt:T];%Creating a mesh grid on the time domain of step dt
N = size(t,2)
u = zeros(1,N);%vector that will store the approximate solution
e = zeros(1,N);%vector that will store the exact solution
u(1) = I;
u(2) = u(1) - 0.5*dt*dt*w*w*u(1);
for i = 2:(N-1) %loop to calculate approximate colution
    u(i+1) = 2*u(i) - u(i-1) - dt*dt*w*w*u(i);
end
for j = 1:N %loop to compute exact solution and store it in vector e
    x = t(j);
    e(j) = I*cos(w*x)
end
error_mat = e-u
maxnorm = max(error_mat)
l2norm = sqrt(dt*sum(error_mat.^2))
f1 = figure
plot(t,u,'b')%plots approximate solution
hold on
plot(t,e,'r' )%plot exact solution
legend('Approximation', 'exact')
xlabel('t')
ylabel('u')
hold off
saveas(f1,'verlet_FD_graph.png')
end

