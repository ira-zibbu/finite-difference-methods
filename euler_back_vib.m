function [u,e,error_mat,t] = euler_back_vib(I,w,dt,T)
%%Solves the equation u_tt = -w^2 u, the standard vibration ODE using the
%%euler backward method 
t = [0:dt:T];%Creating a mesh grid on the time domain of step dt
N = size(t,2)
u = zeros(1,N);%vector that will store the approximate solution
e = zeros(1,N);%vector that will store the exact solution
u(1) = I;%initial values
u(2) = 2*u(1)/(2+dt*dt*w*w);
for i = 2:N-1
    u(i+1) = (2*u(i)-u(i-1))/(1 + w*w*dt*dt);
end
for j = 1:N %loop to compute exact solution and store it in vector e
    x = t(j);
    e(j) = I*cos(w*x);
end
error_mat = e-u;
maxnorm = max(error_mat)
l2norm = sqrt(dt*sum(error_mat.^2))
f3 = figure
plot(t,u,'b')%plots approximate solution
hold on
plot(t,e,'r' )%plot exact solution
legend('Approximation', 'exact')
xlabel('t')
ylabel('u')
hold off
saveas(f3,'vib_euler_back.png')
end
