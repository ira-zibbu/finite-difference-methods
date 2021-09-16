function [U,E,error_mat,x,t] = diffusion_euler_fwd_sine(dt,dx,a,kp,Q,T)
%%The diffusion equation u_t = a*u_xx  where a is a scalar, x is in the closed interval [0,pi] with initial 
%%conditions u(x,0) = I(x), u(0,t) = 0, u(1,t) = 0. k,Q are scalars and can
%%be any number 

x = [0:dx:pi]; %creating mesh in space domain of step size dx
t = [0:dt:T]; %creating mesh in the time domain of step size dt
s1 = size(x);
s2 = size(t); 
Nx = s1(1,2); %number of elements in x
Nt = s2(1,2); %number of elements in t
u1 = zeros(1,Nx); %creating a blank row vector to store the solution for a time level n
u2 = zeros(1,Nx); %creating a blank row vector to store the solution for a time level n+1
U = zeros(Nt,Nx); %to store the total solution
F = a*dt/(dx^2);
b = a*kp*kp;
if F >= 0.5 %stability criteria
    error('Error F value not small enough');
end
syms z;
I(z) = Q*sin(kp*z); %initial function
for i = 1:Nx %loop to assign the initial value to our solution at time level n =0
    zi = dx*(i-1);
    u1(i)= I(zi);
end
U(1,:) = u1; %assiging the solution of n= 0 time level to our solution matrix
for j = 1:(Nt-1) %loop that passes over all time levels
    for k  = 2:(Nx-1) %loop that passes over all space points for a given time level
        u2(k) = u1(k) + F*(u1(k+1) - 2*u1(k) + u1(k-1)); %FD formulation of solution
    end
    u2(1) = 0; %inital condition u(0,t) = 0
    u2(Nx) = 0;%intial condition u(1,t) = 0
    U((j+1),:)= u2;
    u1 = u2;
end
E = zeros(Nt,Nx); %empty matrix to store the exact solution
for i = 1:Nt
    ti = dt*(i-1);
    for j = 1:Nx
        xi = dx*(j-1);
        E(i,j) = Q*exp(-b*ti)*sin(kp*xi);
    end
end
error_mat = E - U; %matrix that computes the error
[A, B] = meshgrid(x,t);
fig_approx = figure;%plots approximate solution
surf(A,B,U) 
title ('Approximate')
xlabel ('x')
ylabel ('t')
saveas(fig_approx, 'approx_euler_fwd_diff.png');
fig_exact = figure;%plots eact solution
surf(A,B,E)
title ('exact');
xlabel ('x')
ylabel ('t')
saveas(fig_exact, 'exact_euler_fwd_diff.png');
end


    