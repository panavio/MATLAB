function [x,y]=DE2_krisanti(t0,tN,y0,y1,h, p,q,g) 
x=t0:h:tN;

y=zeros(1,length(x)); % y solution
fy=zeros(1, length(x)); %first order y'
y(1)=y0;
fy(1)=y1;

for i=2:length(x)
    y(i)=y(i-1)+h*fy(i-1); %euler's method for first order ODE
    sec=g(x(i-1))-p(x(i-1))* fy(i-1) - q(x(i-1)) * y(i-1);
    fy(i)=fy(i-1)+h.*sec; %euler's method for second order ODE
end

end
