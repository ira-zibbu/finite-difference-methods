function [u,v] = KS_explicit (m,c,dx,dt,T)
%solve the minimal keller-segel model. m is parameter \mu, c is \chi, dx is
%the space mesh size, dt is the mesh size for time and T is the final
%time. The initial conditions are u(x,o) = (x^2)/2 - (x^3)/3 and u_x(x,t) =
%v_x(x,0) = 0 for x = 1 and 0.
t = [0:dt:T]; %creating mesh points in time domain
x = [0:dx:1];%creating mesh points in the space domain
Nx = size(x,2); %no of mesh points in time
Nt = size(t,2);% no of mesh points in space
u = zeros(Nt,Nx); %stores solution for cell concentrations
v = zeros(Nt,Nx); %stores solution for attractnt concentration
F = dt/(dx^2);
for i = 1:Nx %loop to set the initial conditions at time level 0 
    u(1,i) = 1;
    v(1,i) = 1 + 0.1*exp(-10*(x(i)^2));
end
for j = 2:Nt %loop to compute the approximate solution 
    for k = 1:Nx
    if k == 1 
        u(j,k)=m*F*(2*u(j-1,k+1) - u(j-1,k)) - c*F*(1/dx)*(u(j-1,k+1) -u(j-1,k))*(2*v(j-1,k+1) - 2*v(j-1,k)) + u(j-1,k);
        v(j,k) = 2*F*v(j-1,k+1)+ (dt-1-2*F)*v(j-1,k) + dt*u(j-1,k);
    end
    if k == Nx
        u(j,k)=m*F*(2*u(j-1,k-1) - u(j-1,k)) - c*F*(1/dx)*(u(j-1,k-1) -u(j-1,k))*(2*v(j-1,k-1) - 2*v(j-1,k)) + u(j-1,k);
        v(j,k) = 2*F*v(j-1,k-1)+ (dt-1-2*F)*v(j-1,k) + dt*u(j-1,k);
    end
    if ~(k == 1 | k == Nx)
         u(j,k)=m*F*(u(j-1,k+1) - u(j-1,k) + u(j-1,k-1)) - c*F*(1/dx)*(u(j-1,k+1) -u(j-1,k))*(v(j-1,k+1) - 2*v(j-1,k) + v(j-1,k-1)) + u(j-1,k);
        v(j,k) = F*v(j-1,k+1)+ (dt-1-2*F)*v(j-1,k) + dt*u(j-1,k) + F*v(j-1,k-1);
    end
    end
end
[A, B] = meshgrid(x,t);
fig_u = figure;%plots approximate solution
surf(A,B,u) 
title ('cell concentration')
xlabel ('x')
ylabel ('t')
saveas(fig_u,'KS_explicit_u.png');
fig_v = figure;%plots approximate solution
surf(A,B,v) 
title ('attractant concentration')
xlabel ('x')
ylabel ('t')
saveas(fig_v,'KS_explicit_v.png');
fig_ut = figure;
plot(x,u(1,:),'r')
hold on
title ('attractant concentration at different time level')
xlabel ('t')
ylabel ('u')
plot(x,u(5,:),'b')
plot(x,u(10,:),'g')
legend('1st time level', '2nd time level','3rd time level')
hold off
sum_u = sum(u') %stores the total sum of cell mass at each time level
sum_v = sum(v')
end
        