clear; clc; close all;

%% Given data
th0 = -5;      % initial angle [deg]
thf = 80;      % final angle   [deg]
tf  = 4;       % total time    [s]
D   = thf - th0;   % total displacement = 85 deg

t = linspace(0, tf, 1000);

%% Linear trajectory with parabolic blends 
% Free design choice: blend time tb (must satisfy 0 < tb <= tf/2).
tb  = 1;                      % blend time [s]
w   = D/(tf - tb);            % constant cruise velocity [deg/s]
alp = w/tb;                   % blend acceleration       [deg/s^2]

fprintf('--- LSPB parameters ---\n');
fprintf('Blend time          tb    = %.4f s\n', tb);
fprintf('Cruise velocity     omega = %.4f deg/s\n', w);
fprintf('Blend acceleration  alpha = %.4f deg/s^2\n', alp);
fprintf('Position at end of 1st blend = %.4f deg\n', th0 + 0.5*alp*tb^2);

th_l = zeros(size(t)); w_l = zeros(size(t)); a_l = zeros(size(t));
for k = 1:numel(t)
    tk = t(k);
    if tk <= tb                          % 1st parabolic blend
        th_l(k) = th0 + 0.5*alp*tk^2;
        w_l(k)  = alp*tk;
        a_l(k)  = alp;
    elseif tk <= tf - tb                 % linear segment
        th_l(k) = th0 + 0.5*alp*tb^2 + w*(tk - tb);
        w_l(k)  = w;
        a_l(k)  = 0;
    else                                 % 2nd parabolic blend
        th_l(k) = thf - 0.5*alp*(tf - tk)^2;
        w_l(k)  = alp*(tf - tk);
        a_l(k)  = -alp;
    end
end

figure('Name','Linear with parabolic blends');
subplot(3,1,1); plot(t, th_l, 'b', 'LineWidth', 1.5); grid on;
ylabel('\theta [deg]'); title('Linear trajectory with parabolic blends');
subplot(3,1,2); plot(t, w_l, 'r', 'LineWidth', 1.5); grid on;
ylabel('\omega [deg/s]');
subplot(3,1,3); plot(t, a_l, 'k', 'LineWidth', 1.5); grid on;
ylabel('\alpha [deg/s^2]'); xlabel('Time [s]');

%% Cubic trajectory
% theta(t) = a0 + a1 t + a2 t^2 + a3 t^3, with zero start/end velocity
a0 = th0;
a1 = 0;
a2 =  3*D/tf^2;
a3 = -2*D/tf^3;

fprintf('\n--- Cubic coefficients ---\n');
fprintf('a0 = %.4f\na1 = %.4f\na2 = %.4f\na3 = %.5f\n', a0, a1, a2, a3);

th_c = a0 + a1*t + a2*t.^2 + a3*t.^3;
w_c  = a1 + 2*a2*t + 3*a3*t.^2;
a_c  = 2*a2 + 6*a3*t;

figure('Name','Cubic trajectory');
subplot(3,1,1); plot(t, th_c, 'b', 'LineWidth', 1.5); grid on;
ylabel('\theta [deg]'); title('Cubic trajectory');
subplot(3,1,2); plot(t, w_c, 'r', 'LineWidth', 1.5); grid on;
ylabel('\omega [deg/s]');
subplot(3,1,3); plot(t, a_c, 'k', 'LineWidth', 1.5); grid on;
ylabel('\alpha [deg/s^2]'); xlabel('Time [s]');


