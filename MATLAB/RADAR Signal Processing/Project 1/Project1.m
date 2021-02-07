%%%% John Claus %%%%
%%%% ECE538     %%%%
%%%% Project 1  %%%%
%%%% 11/08/2019 %%%%

clear all
close all

%%% Problem 1 %%%
% Initial values
G_db = 20;      %db
RCS_db = 30;    %dbsm
Noise_db = 3;   %db
Lsys_db = 0;    %db   
Latm_db = 0;    %db
Pt_db = 10;     %dbW
R = 10000;      %m
c = physconst('LightSpeed');    %mps
Fc = 1e9;       %Hz
Fs = 4*Fc;      %Hz
PW = 1e-6;      %m

% Convert inital values from db
G = 10^(G_db/10);
RCS = 10^(RCS_db/10);
Lsys = 10^(Lsys_db/10);
Latm = 10^(Latm_db/10);
Pt = 10^(Pt_db/10);

% Calculate constants
lambda = c / Fc;
T = 4 * PW;
Ts = 1/Fs;
N = round(T*Fs);
t = (0:(N-1)).*Ts;
LOC = Ts*(round(N/2));

% Calculate Barker Code-13 constants
% +++++--++-+-+
Chirp = PW / 13;
Chirp_L = Chirp / Ts;
Chirp_N = round(Chirp / Ts);

% Calculate Received Power
Pr = (Pt*(G^2)*(lambda^2)*RCS)/(((4*pi)^3)*(R^4)*Lsys*Latm);
sigma_sq = Pr*(10^(Noise_db/10));

% Create output waveforms
x_out = [exp(1j.*2.*pi.*Fc*t(1:round(N/4))) zeros(1,round((3*N)/4))]; 
x_out_BC = [exp(1j.*2.*pi.*Fc*t(1:5*Chirp_N)) ...
    exp(1j.*pi.*Fc*t((5*Chirp_N)+1:7*Chirp_N)) ...
    exp(1j.*2.*pi.*Fc*t((7*Chirp_N)+1:9*Chirp_N)) ...
    exp(1j.*pi.*Fc*t((9*Chirp_N)+1:10*Chirp_N)) ...
    exp(1j.*2.*pi.*Fc*t((10*Chirp_N)+1:11*Chirp_N)) ...
    exp(1j.*pi.*Fc*t((11*Chirp_N)+1:12*Chirp_N)) ...
    exp(1j.*2.*pi.*Fc*t((12*Chirp_N)+1:13*Chirp_N)) ...
    zeros(1,N-(13*Chirp_N))];
    
% Create input waveforms
x = sqrt(Pr)*[zeros(1,round(N/2)) exp(1j.*2.*pi.*Fc*t(1:round(N/4))) ...
    zeros(1,round(N/4))]; 
x_BC = sqrt(Pr)*[zeros(1,round(N/2)) ...
    exp(1j.*2.*pi.*Fc*t(1:5*Chirp_N)) ...
    exp(1j.*pi.*Fc*t((5*Chirp_N)+1:7*Chirp_N)) ...
    exp(1j.*2.*pi.*Fc*t((7*Chirp_N)+1:9*Chirp_N)) ...
    exp(1j.*pi.*Fc*t((9*Chirp_N)+1:10*Chirp_N)) ...
    exp(1j.*2.*pi.*Fc*t((10*Chirp_N)+1:11*Chirp_N)) ...
    exp(1j.*pi.*Fc*t((11*Chirp_N)+1:12*Chirp_N)) ...
    exp(1j.*2.*pi.*Fc*t((12*Chirp_N)+1:13*Chirp_N)) ...
    zeros(1,N -((13*Chirp_N)+round(N/2)))];

% Create noise for received waveforms
x_noise = sqrt(0.5.*sigma_sq).*(randn(size(x))+1j.*randn(size(x)));
x_noise_BC = sqrt(0.5.*sigma_sq).*(randn(size(x_BC))+1j.*randn(size(x_BC)));

y = x + x_noise;
y_BC = x_BC + x_noise_BC;

% Matched Filtering - No modulation
xcorOut = xcorr(x_out,y);
numberOfLags = length(xcorOut);
lagVector = ((1:(2*N-1)));
xcorTimeVec = lagVector.*Ts;

% Matched Filtering - No modulation - No Noise
xcorOut_NN = xcorr(x_out,x);
numberOfLags_NN = length(xcorOut_NN);
lagVector_NN = ((1:(2*N-1)));
xcorTimeVec_NN = lagVector_NN.*Ts;

% Matched Filtering - Barker Code 13
xcorOut_BC = xcorr(x_out_BC,y_BC);
numberOfLags_BC = length(xcorOut_BC);
lagVector_BC = ((1:(2*N-1)));
xcorTimeVec_BC = lagVector_BC.*Ts;

% Matched Filtering - Barker Code 13 - No Noise
xcorOut_BC_NN = xcorr(x_out_BC,x_BC);
numberOfLags_BC_NN = length(xcorOut_BC_NN);
lagVector_BC_NN = ((1:(2*N-1)));
xcorTimeVec_BC_NN = lagVector_BC_NN.*Ts;

% Matched Filtering - Noise Only 
xcorOut_NO = xcorr(x_out,x_noise);
numberOfLags_NO = length(xcorOut_NO);
lagVector_NO = ((1:(2*N-1)));
xcorTimeVec_NO = lagVector_NO.*Ts;

% Plot Unmodulated Signal and Matched Filter Output
figure
subplot(2,1,1)
plot(t./1e-6,y);
ylim([0 3e-6])
ylabel('Amplitude')
xlabel('Time Samples (\musec)')
grid on;title('Receive Signal')
hold on
plot([LOC, LOC]./1e-6, [0 max(abs(y))],'r')
subplot(2,1,2)
plot(xcorTimeVec./1e-6,abs(xcorOut))
ylabel('Cross Correlation Amplitude')
xlabel('Delay Time (\musec)')
grid on;
title('Matched Filter Output')
hold on;
plot([LOC LOC]./1e-6,[0 max(abs(xcorOut))],'r')
hold off;

% Plot Unmodulated Signal and Matched Filter Output - No Noise
figure
subplot(2,1,1)
plot(t./1e-6,x);
ylim([0 3e-6])
ylabel('Amplitude')
xlabel('Time Samples (\musec)')
grid on;title('Receive Signal')
hold on
plot([LOC, LOC]./1e-6, [0 max(abs(x))],'r')
subplot(2,1,2)
plot(xcorTimeVec_NN./1e-6,abs(xcorOut_NN))
ylabel('Cross Correlation Amplitude')
xlabel('Delay Time (\musec)')
grid on;
title('Matched Filter Output - No Noise')
hold on;
plot([LOC LOC]./1e-6,[0 max(abs(xcorOut_NN))],'r')
hold off;

% Plot Barker Code - 13 Signal and Matched Filter Output
figure
subplot(2,1,1)
plot(t./1e-6,y_BC);
ylim([0 3e-6])
ylabel('Amplitude')
xlabel('Time Samples (\musec)')
grid on;title('Receive Signal BC-13')
hold on
plot([LOC, LOC]./1e-6, [0 max(abs(y_BC))],'r')
subplot(2,1,2)
plot(xcorTimeVec_BC./1e-6,abs(xcorOut_BC))
ylabel('Cross Correlation Amplitude')
xlabel('Delay Time (\musec)')
grid on;
title('Matched Filter Output BC-13')
hold on;
plot([LOC LOC]./1e-6,[0 max(abs(xcorOut_BC))],'r')
hold off;

% Plot Barker Code - 13 Signal and Matched Filter Output - No Noise
figure
subplot(2,1,1)
plot(t./1e-6,x_BC);
ylim([0 3e-6])
ylabel('Amplitude')
xlabel('Time Samples (\musec)')
grid on;title('Receive Signal BC-13')
hold on
plot([LOC, LOC]./1e-6, [0 max(abs(x_BC))],'r')
subplot(2,1,2)
plot(xcorTimeVec_BC_NN./1e-6,abs(xcorOut_BC_NN))
ylabel('Cross Correlation Amplitude')
xlabel('Delay Time (\musec)')
grid on;
title('Matched Filter Output BC-13')
hold on;
plot([LOC LOC]./1e-6,[0 max(abs(xcorOut_BC_NN))],'r')
hold off;

% Plot Noise Only and Matched Filter Output 
figure
subplot(2,1,1)
plot(t./1e-6,x_noise);
ylim([0 3e-6])
ylabel('Amplitude')
xlabel('Time Samples (\musec)')
grid on;title('Receive Signal Noise Only')
hold on
plot([LOC, LOC]./1e-6, [0 max(abs(x_noise))],'r')
subplot(2,1,2)
plot(xcorTimeVec_NO./1e-6,abs(xcorOut_NO))
ylabel('Cross Correlation Amplitude')
xlabel('Delay Time (\musec)')
grid on;
title('Matched Filter Output Noise Only')
hold on;
plot([LOC LOC]./1e-6,[0 max(abs(xcorOut_NO))],'r')
hold off;

%%% Problem 2 %%%
% Define Initial Variables and Zero Lists
max_iter = 100000;
match_sig_op = zeros(1,max_iter);
match_noise_op = zeros(1,max_iter);
Pfa = zeros(1,20);
Pd = zeros(1,20);
thr_graph = zeros(1,20);

% Creates List of Values at Expected 0th Matched Filter Output
for iter = 1:max_iter
    x = sqrt(Pr)*[zeros(1,round(N/2)) exp(1j.*2.*pi.*Fc*t(1:round(N/4))) ...
        zeros(1,round(N/4))]; 
    x_noise = sqrt(0.5.*sigma_sq).*(randn(size(x))+1j.*randn(size(x)));
    y = x + x_noise;
    xcorOut = xcorr(x_out,y);
    xcorOut_NO = xcorr(x_out,x_noise);
    index = round(length(xcorOut)/4);
    noise_op = abs(xcorOut_NO(index));
    signal_op = abs(xcorOut(index));
    match_sig_op(1,iter) = signal_op;
    match_noise_op(1,iter) = noise_op;
end

% Calculates Pd and Pfa from the Matched Filter Output Data
iter2 = 1;
for THR = 1:5:100
    noisePwr = var(real(match_noise_op)).*2;
    threshold = noisePwr*THR;
    Pfa(iter2) = sum(abs(match_noise_op).^2>threshold)/max_iter;
    Pd(iter2)  = sum(abs(match_sig_op).^2>threshold)/max_iter;
    thr_graph(iter2) = THR;
    iter2 = iter2 + 1;
end

% Plots Pd vs Pfa data
figure
plot(thr_graph,Pfa,thr_graph,Pd);
legend('Pfa','Pd')
ylim([0 1.5])
ylabel('Amplitude')
xlabel('THR value')
grid on;title('Pfa vs. Pd')





