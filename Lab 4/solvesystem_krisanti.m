function [x,y]=solvesystem_krisanti(f,g, t0, tN, x0, h) %f,start, endpoint, initial condition, stepsize

x=t0:h:tN;
%disp(x);
size(x)
y=zeros(2, length(x));
y(1,1)=x0(1);
y(2,1)=x0(2);
num=size(x);
for i=1:(num(2)-1)
    kf=f(x(i),y(1,i),y(2,i)); %y(tn+1) = y(n) + (1/2) * h * (f(tn,yn) + f(tn+h, yn+h*f(tn,yn))
    kg=g(x(i),y(1,i),y(2,i));
    x1=y(1,i)+h.*kf;
    x2=y(2,i)+h.*kg;
    y(1, i+1)=(kf+f(x(i+1),x1,x2)).*(1/2).*h   +  y(1,i);
    y(2, i+1)=(kg+g(x(i+1),x1,x2)).*(1/2).*h   +  y(2,i);
end
[ox,oy]=ode45(@(t,y) [f(t, y(1), y(2)); g(t, y(1), y(2))], [t0,tN] ,x0);
%plot(x,y)

%legend("mine", "ode45")
end

%x1' = x1/2 - 2*x2, x2' = 5*x1 - x2
%with initial condition x(0)=(1,1).
%Use your method from Exercise 1 to approximate the solution from t=0 to t=4*pi with step size h=0.05.
 
