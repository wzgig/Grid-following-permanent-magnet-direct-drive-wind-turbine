%% test_ode45_ok.m
% 测试 ode45 是否正常工作的脚本
% 1) 指数衰减 y'=-y (解析解 y=e^{-t})
% 2) 简谐振子 x''+x=0 (能量应近似守恒)
% 3) 事件检测：y 下降到 0.5 时停止

clear; clc; close all;

fprintf('==== ODE45 Sanity Check ====\n');

%% Test 1: y' = -y, y(0)=1  (exact: y=exp(-t))
f1 = @(t,y) -y;
tspan1 = [0 5];
y0_1 = 1;

% 默认容差
optsA = odeset('RelTol',1e-6,'AbsTol',1e-9);
[tA,yA] = ode45(f1,tspan1,y0_1,optsA);
yExactA = exp(-tA);
errInfA = max(abs(yA - yExactA));

% 更严格容差（误差应更小）
optsB = odeset('RelTol',1e-9,'AbsTol',1e-12);
[tB,yB] = ode45(f1,tspan1,y0_1,optsB);
yExactB = exp(-tB);
errInfB = max(abs(yB - yExactB));

fprintf('\n[Test 1] y''=-y, y(0)=1, t in [0,5]\n');
fprintf('  Max abs error (RelTol=1e-6, AbsTol=1e-9):  %.3e\n', errInfA);
fprintf('  Max abs error (RelTol=1e-9, AbsTol=1e-12): %.3e\n', errInfB);
fprintf('  Error ratio (A/B): %.2f (通常应 > 1)\n', errInfA/errInfB);

figure('Name','Test1: Exponential Decay');
plot(tA,yA,'o-', tA,yExactA,'-','LineWidth',1.2);
grid on; xlabel('t'); ylabel('y');
legend('ode45','exact e^{-t}','Location','best');
title('Test 1: y''=-y (ode45 vs exact)');

%% Test 2: Simple harmonic oscillator
% x'' + x = 0  =>  y1=x, y2=x'
% y1' = y2
% y2' = -y1
f2 = @(t,y) [y(2); -y(1)];
tspan2 = [0 20*pi];
y0_2 = [1; 0]; % x(0)=1, v(0)=0 => exact x=cos(t)

opts2 = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t2,y2] = ode45(f2,tspan2,y0_2,opts2);

x = y2(:,1); v = y2(:,2);
E = 0.5*(v.^2 + x.^2);          % 能量（理论上恒定 0.5）
E0 = E(1);
relDrift = max(abs(E - E0))/E0; % 相对漂移

fprintf('\n[Test 2] SHO x''''+x=0, t in [0, 20*pi]\n');
fprintf('  Energy relative drift: %.3e (应当很小)\n', relDrift);

figure('Name','Test2: SHO');
subplot(2,1,1);
plot(t2,x,'LineWidth',1.2); grid on;
xlabel('t'); ylabel('x'); title('x(t) from ode45');

subplot(2,1,2);
plot(t2,E,'LineWidth',1.2); grid on;
xlabel('t'); ylabel('Energy'); title('Energy drift (should be near constant)');

%% Test 3: Event detection (stop when y hits 0.5)
optsE = odeset(optsA,'Events',@hitHalfEvent);

[tE,yE,tEvent,yEvent,iEvent] = ode45(f1,[0 10],1,optsE);

tTrue = log(2);
tErr = abs(tEvent - tTrue);

fprintf('\n[Test 3] Event: stop when y(t)=0.5\n');
fprintf('  Detected tEvent: %.12f, true ln(2): %.12f, abs error: %.3e\n', ...
    tEvent, tTrue, tErr);

figure('Name','Test3: Event');
plot(tE,yE,'o-','LineWidth',1.2); hold on;
yline(0.5,'--');
xline(tEvent,'--');
grid on; xlabel('t'); ylabel('y');
title('Event detection: stop at y=0.5');

% ---- local event function (place at end of file) ----
function [value, isterminal, direction] = hitHalfEvent(t,y)
    value = y - 0.5;    % 触发条件 value=0
    isterminal = 1;     % 触发后停止积分
    direction  = -1;    % 只检测从上往下穿过(下降)
end

%% Simple PASS/FAIL rule of thumb (你可以按需要调整阈值)
pass1 = (errInfB < 1e-6) && (errInfA/errInfB > 1.5);
pass2 = (relDrift < 1e-5);
pass3 = (tErr < 1e-6);

fprintf('\n==== Summary ====\n');
fprintf('  Test1 (accuracy & tighter tol improves): %s\n', tfstr(pass1));
fprintf('  Test2 (energy drift small):             %s\n', tfstr(pass2));
fprintf('  Test3 (event time accurate):            %s\n', tfstr(pass3));

if pass1 && pass2 && pass3
    fprintf('==> ode45 看起来工作正常。\n');
else
    fprintf('==> 有测试未通过：可能是设置/版本/函数写法问题，建议检查。\n');
end

%% helper
function s = tfstr(tf)
    if tf, s = 'PASS'; else, s = 'FAIL'; end
end
