% EXTRACTING DATA FROM TABLE
t = scope112(:,1); % Time data
u = scope112(:,2); % Input signal data
y1 = scope112(:,3); % Output signal without zero
y2 = scope112(:,4); % Output signal with zero

% DISPLAYING DATA
plot(t, [u, y1+3, y2+6]), shg
title('Display of acquired data')
lgd = legend('Input: u', 'System output without zero: y1', 'System output with zero: y2');
title(lgd, 'Legend')
xlabel('Time [s]');
ylabel('Amplitude [V]');

%% NON-PARAMETRIC IDENTIFICATION
figure;
plot(t, [u, y1])
title('System without zero')
lgd = legend('Input: u', 'Output1: y1');
title(lgd, 'Legend')
xlabel('Time [s]');
ylabel('Amplitude [V]');

% Resonance identification (DO NOT MODIFY!)
u_max = 448; 
y1_max = 450;
u_min = 455;
y1_min = 457;

% Identifying the proportionality factor
k = mean(y1) / mean(u);

% Amplitude of the module at resonance
Mr = ((y1(y1_max) - y1(y1_min)) / (u(u_max) - u(u_min)));

% Identifying damping factor
zeta = sqrt((Mr - sqrt(Mr^2 - 1)) / (2 * Mr));

% Identifying the natural frequency
wr = (pi) / (t(y1_min) - t(y1_max)); 
wn = wr / (sqrt(1 - 2 * (zeta^2)));

% Transfer function obtained
NUM = k * wn^2;
DEN = [1, 2 * zeta * wn, wn^2];
Hs = tf(NUM, DEN);
ysim = lsim(Hs, u, t);

figure;
plot(t, ysim, t, y1);
title('Overlay of acquired signal and non-parametric identified signal')
lgd = legend('Identified signal', 'Acquired signal');
title(lgd, 'Legend')
xlabel('Time [s]');
ylabel('Amplitude [V]');

%% State-space representation
A = [0 1; -wn^2 -2 * zeta * wn];
B = [0; k * wn^2];
C = [1 0];
D = [0];

sys = ss(A, B, C, D);
figure;
ysim = lsim(sys, u, t, [y1(1), (y1(2) - y1(1)) / (t(2) - t(1))]);
plot(t, ysim, t, y1), shg
title('Display of identified data')
lgd = legend('Identified y', 'Acquired y');
title(lgd, 'Legend')
xlabel('Time [s]');
ylabel('Amplitude [V]');

% Mean squared error
J = norm(y1 - ysim) / sqrt(length(y1));

% Normalized mean squared error
eMPN = norm(y1 - ysim) / norm(y1 - mean(y1)) * 100;
fprintf("Mean Squared Error J: %f\n", J);
fprintf("Normalized Mean Squared Error eMPN: %f%%\n", eMPN);


%% PARAMETRIC IDENTIFICATION
dt = t(2) - t(1);
data_y1 = iddata(y1, u, dt);
M_y1_armax = armax(data_y1, [2, 2, 2, 0]); % nB = number of zeros +1 
figure;
resid(M_y1_armax, data_y1, 5);
figure;
compare(data_y1, M_y1_armax);
H_y1_armax_d = tf(M_y1_armax.B, M_y1_armax.A, dt, 'variable', 'z^-1');
H_y1_armax_c = minreal(zpk(d2c(H_y1_armax_d, 'zoh')));

%% IV4 Method
M_y1_iv4 = iv4(data_y1, [2, 2, 0]);
figure;
resid(M_y1_iv4, data_y1, 5);
figure;
compare(data_y1, M_y1_iv4);
H_y1_iv4_d = tf(M_y1_iv4.B, M_y1_iv4.A, dt, 'variable', 'z^-1');
H_y1_iv4_c = minreal(zpk(d2c(H_y1_iv4_d, 'zoh')));

%% PARAMETRIC IDENTIFICATION FOR OUTPUT y2
data_y2 = iddata(y2, u, dt);
M_y2_armax = armax(data_y2, [2, 2, 2, 0]);
figure;
resid(M_y2_armax, data_y2, 5);
figure;
compare(data_y2, M_y2_armax);
H_y2_armax_d = tf(M_y2_armax.B, M_y2_armax.A, dt, 'variable', 'z^-1');
H_y2_armax_c = minreal(zpk(d2c(H_y2_armax_d, 'zoh')));

%% IV4 Method for y2
M_y2_iv4 = iv4(data_y2, [2, 2, 0]);
M_y2_iv4 = pem(data_y2, M_y2_iv4);
figure;
resid(M_y2_iv4, data_y2, 5);
figure;
compare(data_y2, M_y2_iv4);
H_y2_iv4_d = tf(M_y2_iv4.B, M_y2_iv4.A, dt, 'variable', 'z^-1');
H_y2_iv4_c = minreal(zpk(d2c(H_y2_iv4_d, 'zoh')));

%% FREQUENCY RESPONSE ANALYSIS
% Identification of modulus

u1_max = 113;
y1_max = 115;
u1_min = 139;
y1_min = 141;
A_o1 = (y1(y1_max)-y1(y1_min))/2;
A_i1 = (u(u1_max)-u(u1_min))/2;
M1 = A_o1/A_i1;
dt1 = t(y1_max)-t(u1_max);
w1 = pi/(t(y1_min)-t(y1_max));

u2_max = 196;
y2_max = 196;
u2_min = 211;
y2_min = 212;
A_o2 = (y1(y2_max)-y1(y2_min))/2;
A_i2 = (u(u2_max)-u(u2_min))/2;
M2 = A_o2/A_i2;
dt2 = t(y2_max)-t(u2_max);
w2 = pi/(t(y2_min)-t(y2_max));

u3_max = 276;
y3_max = 277;
u3_min = 285;
y3_min = 288;
A_o3 = (y1(y3_max)-y1(y3_min))/2;
A_i3 = (u(u3_max)-u(u3_min))/2;
M3 = A_o3/A_i3;
dt3 = t(y3_max)-t(u3_max);
w3 = pi/(t(y3_min)-t(y3_max));

u4_max = 354;
y4_max = 356;
u4_min = 363;
y4_min = 365;
A_o4 = (y1(y4_max)-y1(y4_min))/2;
A_i4 = (u(u4_max)-u(u4_min))/2;
M4 = A_o4/A_i4;
dt4 = t(y4_max)-t(u4_max);
w4 = pi/(t(y4_min)-t(y4_max));


u5_max = 448; 
y5_max = 450;
u5_min = 455;
y5_min = 457;
A_o5 = (y1(y5_max)-y1(y5_min))/2;
A_i5 = (u(u5_max)-u(u5_min))/2;
M5 = A_o5/A_i5;
dt5 = t(y5_max)-t(u5_max);
w5 = pi/(t(y5_min)-t(y5_max));


u6_max = 525; 
y6_max = 528;
u6_min = 531;
y6_min = 534;
A_o6 = (y1(y6_max)-y1(y6_min))/2;
A_i6 = (u(u6_max)-u(u6_min))/2;
M6 = A_o6/A_i6;
dt6 = t(y6_max)-t(u6_max);
w6 = pi/(t(y6_min)-t(y6_max));

u7_max = 725;
y7_max = 728;
u7_min = 729;
y7_min = 732;
A_o7 = (y1(y7_max)-y1(y7_min))/2;
A_i7 = (u(u7_max)-u(u7_min))/2;
M7 = A_o7/A_i7;
dt7 = t(y7_max)-t(u7_max);
w7 = pi/(t(y7_min)-t(y7_max));


%PHASE IDENTIFICATION 
u1_max = 113;
y1_max = 115;
y1_min = 141;
dt1 = t(y1_max)-t(u1_max);
w1_p = pi/(t(y1_min)-t(y1_max));
ph1 = -rad2deg(w1_p*dt1);

u2_max = 194;
y2_max = 196;
y2_min = 212;
dt2 = t(y2_max)-t(u2_max);
w2_p = pi/(t(y2_min)-t(y2_max));
ph2 = -rad2deg(w2_p*dt2);

u3_max = 336;
y3_max = 338;
y3_min = 348;
dt3 = t(y3_max)-t(u3_max);
w3_p = pi/(t(y3_min)-t(y3_max));
ph3 = -rad2deg(w3_p*dt3);

u4_max = 354;
y4_max = 356;
y4_min = 365;
dt4 = t(y4_max)-t(u4_max);
w4_p = pi/(t(y4_min)-t(y4_max));
ph4 = -rad2deg(w4_p*dt4);

u5_max = 488;
y5_max = 491;
y5_min = 497;
dt5 = t(y5_max)-t(u5_max);
w5_p = pi/(t(y5_min)-t(y5_max));
ph5 = -rad2deg(w5_p*dt5);

u6_max = 680;
y6_max = 683;
y6_min = 688;
dt6 = t(y6_max)-t(u6_max);
w6_p = pi/(t(y6_min)-t(y6_max));
ph6 = -rad2deg(w6_p*dt6);

u7_max = 725;
y7_max = 728;
y7_min = 732;
dt7 = t(y7_max)-t(u7_max);
w7_p = pi/(t(y7_min)-t(y7_max));
ph7 = -rad2deg(w7_p*dt7);


 w = logspace(3,4);
[M,Ph] = bode(Hs,w);

figure
subplot(211);
semilogx([w1 w2 w3 w4 w5 w6 w7],20*log10([M1 M2 M3 M4 M5 M6 M7]),'O'); grid;
hold on
semilogx(w, squeeze(20 * log10(M))); grid
title("AMPLITUDE [dB]");
hold on

subplot(212);
semilogx([w1_p w2_p w3_p w4_p w5_p w6_p w7_p],[ph1 ph2 ph3 ph4 ph5 ph6 ph7],'o'); grid;
hold on
semilogx(w, squeeze(Ph)); grid
title("PHASE [deg]");
hold on


%% BODE verification
M6 = 1.2586;
M7 = 0.6207;
w6 = 5.2360e+03;
w7 = 7.8540e+03;
20*log10(M7/M6)*10*(w6/w7);