%% Systems Lab: Systems of ODEs in MATLAB
% In this lab, you will write your own ODE system solver for the Heun method 
% (aka the Improved Euler method), and compare its results to those of |ode45|.
% 
% You will also learn how to save images in MATLAB.
% 
% Opening the m-file lab4.m in the MATLAB editor, step through each part using 
% cell mode to see the results. Compare the output with the PDF, which was generated 
% from this m-file.
% 
% There are four (4) exercises in this lab that are to be handed in on the due 
% date of the lab. Write your solutions in a separate file, including appropriate 
% descriptions in each step. Save the m-files and the pdf-file for Exercise 4 
% and submit them on Quercus.
%% Student Information
%%
% 
%  Student Name: Patricia Nadia Krisanti
%
%%
% 
%  Student Number: 1009669404
%
%% Exercise 1
% Objective: Write your own ODE system solver using the Heun/Improved Euler 
% Method and compare it to |ode45|.
% 
% Details: Consider the system of 2 ODEs:
% 
% |x1'=f(t,x1,x2), x2'=g(t,x1,x2)|
% 
% This m-file should be a function which accepts as variables (t0,tN,x0,h), 
% where t0 and tN are the start and end points of the interval on which to solve 
% the ODE, h is the stepsize, and x0 is a vector for the initial condition of 
% the system of ODEs |x(t0)=x0|. Name the function solvesystem_<UTORid>.m (Substitute 
% your UTORid for <http://utorid UTORid>). You may also want to pass the functions 
% into the ODE the way |ode45| does (check MATLAB labs 2 and 3).
% 
% Your m-file should return a row vector of times and a matrix of approximate 
% solution values (the first row has the approximation for |x1| and the second 
% row has the approximation for |x2|).
% 
% Note: you will need to use a loop to do this exercise. You will also need 
% to recall the Heun/Improved Euler algorithm learned in lectures. 

%{

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
 
%}
%% Exercise 2
% Objective: Compare Heun with an exact solution
% 
% Details: Consider the system of ODEs
% 
% |x1' = x1/2 - 2*x2, x2' = 5*x1 - x2|
% 
% with initial condition |x(0)=(1,1)|.
% 
% Use your method from Exercise 1 to approximate the solution from |t=0| to 
% |t=4*pi| with step size |h=0.05|.
% 
% Compute the exact solution (by hand) and plot both phase portraits on the 
% same figure for comparison.
% 
% Your submission should show the construction of the inline function, the use 
% of your Heun's method to obtain the solution, a construction of the exact solution, 
% and a plot showing both. In the comments, include the exact solution.
% 
% Label your axes and include a legend.

x1=@(t,x1, x2) x1/2-2*x2
x2=@(t,x1,x2) 5*x1 - x2
[x,y]=solvesystem_krisanti(x1, x2, 0, 4*pi, [1,1], 0.05)
plot(y(1,:),y(2,:))
legend('Heun solution to 2 system ODE')
xlabel('x1(t)')
ylabel('x2(t)')
%exact solution for system with initial condition x(0)=(1 1) is 
%c1=1/20
%c2=17/(20sqrt(151))
%s=sqrt(151)
%o=sqrt(151)/4
%x1=c1exp(-t/4).*(3cos(ot)-s*sin(ot)+c2*exp(-t/4)*(s*cos(ot)+3*sin(ot))
%x1=c1exp(-t/4).*(20*cos(ot)+c2*exp(-t/4)*(20*sin(ot))
%% Exercise 3
% Objective: Compare your method with Euler's Method (from |iode|).
% 
% Details: Use |iode| to plot the solution for the same problem with the same 
% step size as on Exercise 2.
% 
% Compare your solution on exercise 2, the exact solution from exercise 2 and 
% the approximation using Euler's method. Plot the solution for Euler's method 
% and make note of any differences.

%euler method using euler.m
fe=@(t,x) [x1(t, x(1), x(2)); x2(t, x(1), x(2))];
t=0:0.05:4*pi
eu=euler(fe, [1;1], t)
plot(eu(1,:), eu(2,:))
%ex(1,:)

%exact solution
c1=1/20
c2=17/(20*sqrt(151))
s=sqrt(151)
o=sqrt(151)/4
ex1=c1.*exp(-t./4).*(3.*cos(o.*t)-s.*sin(o.*t))+c2.*exp(-t./4).*(s.*cos(o.*t)+3.*sin(o.*t))
ex2=c1.*exp(-t./4).*(20.*cos(o.*t))+c2.*exp(-t./4).*(20.*sin(o*t))

%plotting it all together
plot(eu(1,:), eu(2,:),y(1,:),y(2,:),'--', ex1,ex2, 'x')
legend("euler method", "Heun method", "exact solution")
xlabel('x1(t)')
ylabel('x2(t)')


%Both Heun and exact solution shows the 0,0 point. From the initial condition, 
%the exact solution and Heun method approaches 0,0 by moving counterclockwise while the 
%The Euler method did not reach that point, it has a larger error and
%overshoots in the beginning.
%% Saving Images in MATLAB
% To do the following exercises, you will need to know how to output graphics 
% from MATLAB. Create a folder on your Desktop (or elsewhere) to contain the files 
% generated by these exercises. Make this folder the "Current Folder" in the left 
% side of the main MATLAB window. This will ensure that the files output by MATLAB 
% end up in the folder you created.
% 
% To save an image of a phase portrait, use the following steps:
% 
% 1. Get the phase portrait looking the way you want in the |iode| window. 
% 
% 2. Leaving |iode| open, switch to the main MATLAB window.
% 
% 3. Type the command |print -dpng -r300 'filename.png'| in the command window.
% 
% This command will create a PNG graphic called |filename.png| in the current 
% folder. The |-dpng| option tells MATLAB to output the graphic in PNG format; 
% MATLAB also allows output in other formats, such as BMP, EPS, PNG and SVG. The 
% |-r300| option tells MATLAB to set the resolution at 300 dots per inch and can 
% be adjusted if you wish.
%% Exercise 4
% Objective: Analyze phase portraits.
% 
% Details: Compile the results of the following exercises into a single document 
% (e.g. using a word processor) and export it to |PDF| for submission on Quercus. 
% 
% For each of the first-order systems of ODEs 4.1 to 4.10 below, do the following 
% exercises:
% 
% (a) Generate a phase portrait for the system (centre the graph on the equilibrium 
% point at (0,0)). Include a few trajectories.
% 
% (b) Classify the equilibrium on asymptotic stability, and behaviour (sink, 
% source, saddle-point, spiral, center, proper node, improper node) - check table 
% 3.5.1 and figure 3.5.7. Classify also as for clockwise or counterclockwise movement, 
% when relevant.
% 
% (c) Compute the eigenvalues of the matrix (you do not need to show your calculations). 
% Using the eigenvalues you computed, justify part (b).
% 
% To avoid numerical error, you should use Runge-Kutta solver with a step size 
% of |0.05|. Change the display parameters, if necessary, to best understand the 
% phase portrait.
% 
% 4.1. |dx/dt = [2 1; 1 3] x|
% 
% 4.2. |dx/dt = [-2 -1; -1 -3] x|
% 
% 4.3. |dx/dt = [-4 -6; 3 5] x|
% 
% 4.4. |dx/dt = [4 6; -3 -5] x|
% 
% 4.5. |dx/dt = [0 -1; 1 -1] x|
% 
% 4.6. |dx/dt = [0 1; -1 1] x|
% 
% 4.7. |dx/dt = [2 8; -1 -2] x|
% 
% 4.8. |dx/dt = [-2 -8; 1 2] x|
% 
% 4.9. |dx/dt = [-8 5; -13 8] x|
% 
% 4.10. |dx/dt = [8 -5; 13 -8] x|