%Question 6
%Parameter values
c = 0.05; %Fluid flow coefficient, c* (s^-1)
p_infine = 20e6; %Wall rock pressure, p∞ (Pa)
beta = 1e-10; %Compressibility, β (Pa^−1)
v = 1e-6; %Slip rate, v (m s^−1)
d = 10e-6; %Characteristic slip, dc (m)
porosity_ss = 0.03; %Characteristic porosity, φss
porosity_i = 0.025; %Initial porosity, φ0

N=12;
h_value=logspace(-2,1,N); %Generate logarithmically spaced step sizes
pressure_errors=zeros(1,N); %Initialize array to store pressure errors
porosity_errors=zeros(1,N); %Initialize array to store porosity errors


for k=1:N
% discretise
h=h_value(k); % step size
t=0:h:100; % Time vector from 0 to 100 with step size h
S=length(t);

% Initializing an array full of zeros for the solution y. The size is 2 rows and S columns
Y = zeros(2, S);

% Initial conditions
Y(1,1) = 20e6;%initial value of pressure in Pa
Y(2,1) = 0.025;%initial value of porosity

% loop to solve the system of ODEs using the 4th-order Runge-Kutta method
for i=1:S-1
    k1 = f(t(i), Y(:,i), c ,p_infine, beta, v, d, porosity_ss, porosity_i);
    k2 = f(t(i) + h/2, Y(:,i) + k1*h/2, c ,p_infine, beta, v, d, porosity_ss, porosity_i);
    k3 = f(t(i) + h/2, Y(:,i) + k2*h/2, c ,p_infine, beta, v, d, porosity_ss, porosity_i);
    k4 = f(t(i) + h, Y(:,i) + k3*h, c ,p_infine, beta, v, d, porosity_ss, porosity_i);
    Y(:,i+1) = Y(:,i) + h*(k1/6 + k2/3 + k3/3 + k4/6);
end

%Numerical Solutions
y1 = Y(1,:);%Numerical solution of Pressure in Pa
y2 = Y(2,:);%Numerical solution of Porosity

%Exact Solutions
%Exact solution of Pressure in Pa
y1_closed = (((porosity_i - porosity_ss)/beta)/(1- c*d/v)) *(exp(-c*t) - exp(-v*t/d)) + p_infine ;
%Exact solution of Porosity
y2_closed = (porosity_i - porosity_ss)*exp(-v*t/d) + porosity_ss;

% Calculate errors
pressure_error(k)=mean(abs(y1_closed-Y(1,:)));%Pressure error
porosity_error(k)=mean(abs(y2_closed-Y(2,:)));%Porosity error
end

% Plotting
loglog(h_value,pressure_error,'yo','linewidth',2)% Plot pressure error
hold on 
loglog(h_value,porosity_error,'gd','linewidth',2)% Plot porosity error
xlabel('step size')
ylabel('average error')
title('Average Error')
legend('pressure error', 'porosity error')
