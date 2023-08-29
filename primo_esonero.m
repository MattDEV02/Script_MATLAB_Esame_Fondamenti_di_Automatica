syms y(t)

Dy(t) = diff(y, t)
DDy(t) = diff(y, t, 2)
DDDy(t) = diff(y, t, 3)

ode = DDDy(t) + 5 * DDy(t) + 2 * Dy(t) + 6 * y(t) == 1; % 1 per t > 0

cond1 = y(0) == 2;
cond2 = Dy(0) == 0;
cond3 = DDy(0) == 0;

conds = [cond1, cond2, cond3];

ySol(t) = dsolve(ode, conds)

roots = roots(F.Denominator{1})
tau = 1 / abs(roots(2)) 
xt = [0 : 0.01 : 3 * tau]; % t > 0, xt rappresenta il tempo
y(t) = double(ySol(xt));
plot(xt, y(t)), grid

s = tf([1, 0], 1);
F = 1 / (s^3 + 5 * s^2 + 2 * s + 6)

step(F)
stepinfo(F)

rlocus(F)

Kc = 1.1
W = feedback(Kc * F, 1)
step(W)
stepinfo(W)

C = 1 / s
rlocus(C * F)

Kc = 0.5
W = feedback(Kc * C * F, 1)
step(W)
stepinfo(W)

Kp = 3
Kd = 5
Ki = 3
N = 100; %
PID = Kp * ( 1 + Ki / s +Kd * N / ( 1 + N / s))

W2 = feedback(PID * F, 1);
step(W2)
stepinfo(W2)


