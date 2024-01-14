function [t,fin]=f(f, t0, tN, y0, h) %f,start, endpoint, initial condition, stepsize
%calculate euler solution for Z where h=h/2 and there's twice as much
%timesteps as Y
ori=h;
t=[t0]
y=[y0]
z=[y0]
fin=[y0]
tol=1e-2 %this takes foreverrrr (at least on my laptop) so to test increase the tolerance first..
while t(end)<tN
    %t(end)
    %calculate euler solution for Y with 
    k=f(t(end),y(end)); %y(tn+1) = y(n) + (1/2) * h * (f(tn,yn) + f(tn+h, yn+h*f(tn,yn))
    sol=f((t(end)+h)  ,  (y(end)+  h.*k )  );
    y(end+1)=(sol+k).*(1/2).*h   +  y(end);
   

    k=f(t(end),z(end)); %y(tn+1) = y(n) + (1/2) * h * (f(tn,yn) + f(tn+h, yn+h*f(tn,yn))
    sol=f((t(end)+(h/2))  ,  (z(end)+  (h/2).*k )  );
    z(end+1)=(sol+k).*(1/2).*(h/2)   +  z(end);
    

    k=f(t(end),z(end)); %y(tn+1) = y(n) + (1/2) * h * (f(tn,yn) + f(tn+h, yn+h*f(tn,yn))
    sol=f((t(end)+(h/2))  ,  (z(end)+  (h/2).*k )  );
    z(end)=(sol+k).*(1/2).*(h/2)   +  z(end);
    
    
    D=z(end)-y(end);
    if abs(D)<tol
        fin(end+1)=z(end)+D;
        t(end+1)=t(end)+h;
        h=ori;
    elseif abs(D)>=tol
        
        z(end)=[];
        y(end)=[];
        a=tol/abs(D);
        b=min(max(a,0.3),2);
        h = 0.9*h*b;
        
    end

end
    

    s=ode45(f, [t0,tN] ,y0);
    plot(t,z,t,fin,s.x,s.y,'x');
    legend("z","adaptive euler method", "ode45",'Location','NorthWest');
    fin
    
end