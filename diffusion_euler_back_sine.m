function [Total_sol,E,error_mat,x,t] = diffusion_euler_back_sine(dt,dx,a,kp,Q,T)
%%The diffusion equation u_t = a*u_xx  where a is a scalar, x is in the closed interval [0,pi] with initial 
%%conditions u(x,0) = I(x), u(0,t) = 0, u(1,t) = 0. k,Q are scalars and can
%%be any number and are paramters for the exact solution. The euler
%%backward method is an implicit scheme AU=B 

x = [0:dx:pi]; %creating mesh in space domain of step size dx
t = [0:dt:T]; %creating mesh in the time domain of step size dt
s1 = size(x);
s2 = size(t); 
Nx = s1(1,2); %number of elements in x (no. of mesh points)
Nt = s2(1,2); %number of elements in t (no. of mesh points)
Total_sol = zeros(Nt,Nx); %to store the total solution
F = a*dt/(dx^2);
b = a*kp*kp;
B = zeros(1,Nx); %stores the solution to n-1 time level
syms z;
I(z) = Q*sin(kp*z); %initial function
for i = 1:Nx %loop to assign the initial value to our solution at time level n =0
    zi = dx*(i-1);
    B(i)= I(zi);
end
u = zeros(1,Nx);
Total_sol(1,:) = B; %assiging the solution of n= 0 time level to our solution matrix
A = zeros(Nx); %tridiagnol matrix of coefficients
A(1,1) = 1;
A(1,2) = 0;
A(Nx,(Nx-1))=0;
A(Nx,Nx) = 1;
for i = 2:(Nx-1) %loop to out values for A
    A(i,(i-1)) = -F;
    A(i,i) = 2*F+1;
    A(i,(i+1)) = -F;
end
for k = 2:Nt %loop to solve for U at each time level
    U = linsolve(A,B');
    B = U';
    Total_sol(k,:) = U;
end    
%[p,q] = meshgrid(0:0.1:1,0:0.1:1);
%Z = exp(-pi*pi.*q)*sin(pi.*p) + 0.1*exp(-pi*pi*10000.*q)*sin(100*pi.*p);
%surf(p,q,Z)
E = zeros(Nt,Nx); %empty matrix to store the exact solution
for i = 1:Nt
    ti = dt*(i-1);
    for j = 1:Nx
        xi = dx*(j-1);
        E(i,j) = Q*exp(-b*ti)*sin(kp*xi);
    end
end
error_mat = E - Total_sol; %matrix that computer the error
[A, B] = meshgrid(x,t);
fig_approx = figure;%plots approximate solution
surf(A,B,Total_sol) 
title ('Approximate')
xlabel ('x')
ylabel ('t')
saveas(fig_approx,'approx_euler_back.png');
fig_exact = figure;%plots eact solution
surf(A,B,E)
title ('exact');
xlabel ('x')
ylabel ('t')
saveas(fig_exact,'exact_euler_back.png');
end


    