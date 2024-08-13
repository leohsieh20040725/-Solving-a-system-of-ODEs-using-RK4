%Question 4
%Parameter values
c = 0.05; %Fluid flow coefficient, c* (s^-1)
p_infine = 20e6; %Wall rock pressure, p∞ (Pa)
beta = 1e-10; %Compressibility, β (Pa^−1)
v = 1e-6; %Slip rate, v (m s^−1)
d = 10e-6; %Characteristic slip, dc (m)
porosity_ss = 0.03; %Characteristic porosity, φss
porosity_i = 0.025; %Initial porosity, φ0

%Array of times t from 0 to 100 s with a step size h = 1 s
N=101;
tfinal = 100;
t = linspace(0,tfinal,N);
h = t(2) - t(1);

% Initializing an array full of zeros for the solution y. The size is 2 rows and N columns
Y = zeros(2, N);

%Initial conditions for the solution array y
Y(1,1) = 20e6; %initial value of pressure in Pa
Y(2,1) = 0.025;%initial value of porosity


% loop to solve the system of ODEs using the 4th-order Runge-Kutta method
for i=1:N-1
    k1 = f(t(i), Y(:,i), c ,p_infine, beta, v, d, porosity_ss, porosity_i);
    k2 = f(t(i) + h/2, Y(:,i) + k1*h/2, c ,p_infine, beta, v, d, porosity_ss, porosity_i);
    k3 = f(t(i) + h/2, Y(:,i) + k2*h/2, c ,p_infine, beta, v, d, porosity_ss, porosity_i);
    k4 = f(t(i) + h, Y(:,i) + k3*h, c ,p_infine, beta, v, d, porosity_ss, porosity_i);
    Y(:,i+1) = Y(:,i) + h*(k1/6 + k2/3 + k3/3 + k4/6);
end

%Question5
%Numerical Solutions
y1 = Y(1,:);%Numerical solution of Pressure in Pa
y2 = Y(2,:);%Numerical solution of Porosity

%Exact Solutions
%Exact solution of Pressure in Pa
y1_closed = (((porosity_i - porosity_ss)/beta)/(1- c*d/v)) *(exp(-c*t) - exp(-v*t/d)) + p_infine ;
%Exact solution of Porosity
y2_closed = (porosity_i - porosity_ss)*exp(-v*t/d) + porosity_ss;

% Plotting Numerical and Exact Porosity
plot (t, y2, 'ko', 'LineWidth',1);% Numerical Porosity
hold on 
plot (t,y2_closed, 'r-','LineWidth',1);% Exact Porosity
xlabel('Time(s)');
ylabel('Porosity');
legend('Numerical', 'Exact');
title('Porosity vs Time');

% Creating a new figure for Pressure plots
figure()
plot (t, y1, 'ko', 'LineWidth',1);% Numerical Pressure
hold on 
plot (t,y1_closed, 'r-','LineWidth',1);% Exact Pressure
xlabel('Time(s)');
ylabel('Pressure (Pa)');
legend('Numerical', 'Exact');
title('Pressure vs Time');
