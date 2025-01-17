data = readmatrix('scope_2.csv', 'NumHeaderLines', 2);

t = data(:, 1);
u = data(:, 2);
y1 = data(:, 3);
y2 = data(:, 4);

figure
plot(t, [u y1 - 3 y2 - 6])

figure
plot(t, [u y1]); grid
title('Acquired data')
xlabel('Time [s]', 'FontSize', 12)
ylabel('Amplitude [V]', 'FontSize', 12)
legend('u(t)-input', 'y1(t)-output (no zero)', 'FontSize', 10)

u_max = 202;
u_min = 223;
y1_max = 210;
y1_min = 232;

K = mean(y1) / mean(u); % proportionality factor

Mr = ((y1(y1_max) - y1(y1_min)) / (u(u_max) - u(u_min))); % resonance magnitude

T = 2 * (t(y1_min) - t(y1_max)); % period

wr = pi / (t(y1_min) - t(y1_max)); % resonance frequency

zeta = sqrt((Mr - sqrt(Mr^2 - 1)) / 2 / Mr); % damping factor

wn = wr / sqrt(1 - 2 * zeta^2); % natural frequency

NUM = K * wn^2;
DEN = [1, 2 * zeta * wn, wn^2];
H = tf(NUM, DEN);

A = [0, 1; -wn^2, -2 * zeta * wn];
B = [0; K * wn^2];
C = [1, 0];
D = 0;

sys = ss(A, B, C, D);
figure
ysim = lsim(sys, u, t, [y1(1), (y1(2) - y1(1)) / (t(2) - t(1))]);
hold on
plot(t, [y1 ysim]), grid
title('Overlap between measured and simulated signal (2nd order, no zero)', 'FontSize', 11)
xlabel('Time [s]', 'FontSize', 12)
ylabel('Amplitude [V]', 'FontSize', 12)
legend('Output', 'Simulated output', 'FontSize', 10)

J = norm(y1 - ysim) / sqrt(length(y1)); % mean square error = 0.0350
eMPN = norm(y1 - ysim) / norm(y1 - mean(y1)) * 100; % normalized mean square error = 6.1687%

%% Frequency approximation

ph_r = (t(u_max) - t(y1_max)) * wr * 180 / pi;

% Low frequencies

w1 = pi / (t(95) - t(41));
w2 = pi / (t(153) - t(126));
w3 = pi / (t(202) - t(178));

M1 = (y1(47) - y1(101)) / (u(41) - u(95));
M2 = (y1(132) - y1(162)) / (u(126) - u(153));
M3 = (y1(186) - y1(210)) / (u(178) - u(202));

ph1 = (t(41) - t(47)) * w1 * 180 / pi;
ph2 = (t(126) - t(132)) * w2 * 180 / pi;
ph3 = (t(178) - t(186)) * w3 * 180 / pi;

% Medium frequencies

w5 = pi / (t(412) - t(397));
w6 = pi / (t(442) - t(427));
w7 = pi / (t(469) - t(456));

M5 = (y1(408) - y1(423)) / (u(397) - u(412));
M6 = (y1(437) - y1(450)) / (u(427) - u(442));
M7 = (y1(467) - y1(480)) / (u(456) - u(469));

ph5 = (t(398) - t(408)) * w5 * 180 / pi;
ph6 = (t(427) - t(437)) * w6 * 180 / pi;
ph7 = (t(457) - t(467)) * w7 * 180 / pi;

% High frequencies

w8 = pi / (t(765) - t(753));
w9 = pi / (t(787) - t(776));
w10 = pi / (t(809) - t(798));

M8 = (y1(762) - y1(774)) / (u(753) - u(765));
M9 = (y1(785) - y1(795)) / (u(776) - u(787));
M10 = (y1(807) - y1(818)) / (u(798) - u(809));

ph8 = (t(753) - t(762)) * w8 * 180 / pi;
ph9 = (t(776) - t(785)) * w9 * 180 / pi;
ph10 = (t(798) - t(807)) * w10 * 180 / pi;

w = logspace(3, 4);
[M, Ph] = bode(sys, w);

figure
subplot(211)
semilogx([w1 w2 w3 wr w5 w6 w7 w8 w9 w10], 20 * log10([M1 M2 M3 Mr M5 M6 M7 M8 M9 M10]), 'o');
hold on
semilogx(w, squeeze(20 * log10(M))); grid
title('Magnitude diagram')
hold on

subplot(212)
semilogx([w1 w2 w3 wr w5 w6 w7 w8 w9 w10], [ph1 ph2 ph3 ph_r ph5 ph6 ph7 ph8 ph9 ph10], 'o');
hold on
semilogx(w, squeeze(Ph)); grid
title('Phase diagram')
hold on

%% Parametric methods

dt = t(2) - t(1); % sampling period
data_y1 = iddata(y1, u, dt);
data_y2 = iddata(y2, u, dt);

% ARMAX validated by AUTOCORRELATION
armax_y1 = armax(data_y1, [2 1 2 1]);
M_armax_y1 = idpoly(armax_y1);
Hd_armax_y1 = tf(M_armax_y1.B, M_armax_y1.A, dt, 'variable', 'z^-1');
Hd_armax_y1_s = minreal(zpk(Hd_armax_y1));
Hc_armax_y1 = minreal(zpk(d2c(Hd_armax_y1, 'zoh')));

figure
subplot(211)
resid(M_armax_y1, data_y1, 5);
title('Resid ARMAX y1');
subplot(212)
compare(data_y1, M_armax_y1);
title('Compare ARMAX y1');

% ARMAX validated by AUTOCORRELATION
armax_y2 = armax(data_y2, [2 2 2 1]);
M_armax_y2 = idpoly(armax_y2);
Hd_armax_y2 = tf(M_armax_y2.B, M_armax_y2.A, dt, 'variable', 'z^-1');
Hd_armax_y2_s = minreal(zpk(Hd_armax_y2));
Hc_armax_y2 = minreal(zpk(d2c(Hd_armax_y2, 'zoh')));

figure
subplot(211)
resid(M_armax_y2, data_y2, 5);
title('Resid ARMAX y2');
subplot(212)
compare(data_y2, M_armax_y2);
title('Compare ARMAX y2');

% OE validated by CROSS-CORRELATION
oe_y1 = oe(data_y1, [2 2 1]);
M_oe_y1 = idpoly(oe_y1);
Hd_oe_y1 = tf(M_oe_y1.B, M_oe_y1.F, dt, 'variable', 'z^-1');
Hd_oe_y1_s = minreal(zpk(Hd_oe_y1));
Hc_oe_y1 = minreal(zpk(d2c(Hd_oe_y1, 'zoh')));

figure
subplot(211)
resid(M_oe_y1, data_y1, 5);
title('Resid OE y1');
subplot(212)
compare(data_y1, M_oe_y1);
title('Compare OE y1');

% OE validated by CROSS-CORRELATION
oe_y2 = oe(data_y2, [2 3 1]);
M_oe_y2 = idpoly(oe_y2);
Hd_oe_y2 = tf(M_oe_y2.B, M_oe_y2.F, dt, 'variable', 'z^-1');
Hd_oe_y2_s = minreal(zpk(Hd_oe_y2));
Hc_oe_y2 = minreal(zpk(d2c(Hd_oe_y2, 'zoh')));

figure
subplot(211)
resid(M_oe_y2, data_y2, 5);
title('Resid OE y2');
subplot(212)
compare(data_y2, M_oe_y2);
title('Compare OE y2');
