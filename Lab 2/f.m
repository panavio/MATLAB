function [x,y]=f(f, t0, tN, y0, h) %f,start, endpoint, initial condition, stepsize

x=t0:h:tN;
disp(x);
y=zeros(size(x));
y(1)=y0;
num=size(x);
for i=1:(num(2)-1)
    k=f(x(i),y(i)); %y(tn+1) = y(n) + (1/2) * h * (f(tn,yn) + f(tn+h, yn+h*f(tn,yn))
    sol=f((x(i)+h)  ,  (y(i)+  h.*k )  );
    y(i+1)=(sol+k).*(1/2).*h   +  y(i);
end
y(1)=y0;
s=ode45(f, [t0,tN] ,y0);
plot(x,y,s.x,s.y,'x');
legend("mine", "ode45")
end

 
