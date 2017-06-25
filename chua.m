R = 1000;
C1 = 10e-9;
C2 = 100e-9;
L = 18e-3;
xlinha = @(t,x) [1/(R*C1)*(x(2)-x(1))-1/C1*g(x(1)); 1/(R*C2)*(x(1)-x(2))+1/C2*x(3); -1/L*x(2)];
x0 = [-0.5;-0.2;0];
[t,x] = ode45(xlinha , [0 0.05],x0);
plot3(x(:,1),x(:,2),x(:,3))
xlabel('V1')
ylabel('V2')
zlabel('IL')