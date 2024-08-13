%Question 3
% This is a function definition with inputs t, Y, c, p_infine, beta, v, d, porosity_ss, porosity_i
function out = f(t, Y, c ,p_infine, beta, v, d, porosity_ss, porosity_i)



% Extract the pressure (P) and porosity (φ) values from the input vector Y
y1 = Y(1);% Pressure (P) at time t
y2 = Y(2);% Porosity (φ) at time t

% Equation (1): Rate of change of pressure (dp/dt)
f1 = (v/(d*beta))*(y2 - porosity_ss) + c*(p_infine - y1);

% Equation (2): Rate of change of porosity (dφ/dt)
f2 = -(v/d) * (y2 - porosity_ss);

% Combine the results into a column vector
out = [f1;f2];