clear all;
set(0, 'defaultfigurecolor', 'w') % set white ground
T = 1000; 

P = [0.5 0.5 -.3 -.3 -.2 -.2];  % assume primary path
Sec = P;

x = randn(1,T);
y = filter(Sec, 1, x);
S_estx = zeros(1,length(P)); 
S_estw = zeros(1,length(P));     % the weight of Sh(z)
e = zeros(1, T);                 % data buffer for the identification error

% apply LMS algorithm
mu = 0.09;                         
for k = 1 : T                      
    S_estx = [x(k) S_estx(1 : (length(P) - 1))];  
    S_esty = sum(S_estx.*S_estw);	        
    e(k) = y(k) - S_esty;     
    S_estw = S_estw+mu * e(k) * S_estx;   
end

 % check the result
 subplot(2, 1, 1);
 plot([1 : T], e);
 ylabel('Amplitude');
 xlabel('time k');
 legend('error');
 subplot(2,1,2);
 stem(Sec);
 hold on;
 stem(S_estw, '*');
 ylabel('Amplitude');
 xlabel('Numbering of filter tap');
 legend('Coefficients of Sec(z)', 'Coefficients of Sest(z)')

[y,Fs] = audioread('sample.wav');
dt = 1 / Fs;
t = 0 : dt : (length(y) * dt) - dt;
figure;         % plot original voice file in time domain
subplot(2, 2, 1);
plot(t, y);        
xlabel("Time");
ylabel("Amplitude");
title("Original voice file in time domain");

y = y.';
[S, F, T] = stft(y(1,:), Fs, 'Window', kaiser(256,5), 'OverlapLength', 128, 'FFTLength', 256);  % apply stft       
% plot original voice file in frequency domain
subplot(2, 2, 3);
waterfall(F, T, abs(S)');
xlabel("Frequency");
ylabel("Time");
zlabel("Abs()");
title("Original voice file in frequency domain");

% find white noise in voice file
SIZE_S = size(S);
min_S = zeros(SIZE_S(1), SIZE_S(2));
for i = 1 : SIZE_S(1)
    g = 1;
    for k = 1 : 9
        tmp = min(abs(S(i, g : g + int32(SIZE_S(2) / 10) - 1)));
        for j = g : g + int32(SIZE_S(2) / 10) - 1
            min_S(i, j) = tmp;
        end
        g = g + int32(SIZE_S(2) / 10) - 1;
    end
end

% plot white noise in frequency domain
subplot(2, 2, 4);
waterfall(F, T, abs(min_S)');
xlabel("Frequency");
ylabel("Time");
zlabel("Abs()");
title("White noise in frequency domain");

% transfer back to time domain
x = istft(min_S, Fs, 'Window', kaiser(256,5), 'OverlapLength', 128, 'FFTLength', 256);
% plot white noise in time domain
subplot(2, 2, 2);
plot(t,x);
xlabel("Time");
ylabel("Amplitude");
title("White noise in time domain");

X = x.';    % let X beecome white noise in time domain
T = length(X);
P = S_estw;
Yd = filter(P, 1, X);
  
% Initiate the system,
Cx = zeros(1,length(P));       % the state of C(z)
Cw = zeros(1,length(P));       % the weight of C(z)
Sx = zeros(size(Sec));         % the dummy state for the secondary path
err = zeros(1,T);              % data buffer for the control error
Xhx = zeros(1,length(P));      % the state of the filtered x(k)
X_estx = zeros(1,length(P));

% apply the FxLMS algorithm
mu = 100000;                   % step size      
for k = 1:T                        
    Cx = [X(k) Cx(1:length(P) - 1)];           
    Cy = sum(Cx .* Cw);          
    Sx = [Cy Sx(1:length(Sx) - 1)];    
    err(k) = Yd(k) - sum(Sx .* S_estw);  
    X_estx = [X(k) X_estx(1:length(P) - 1)];         
    Xhx = [sum(X_estx .* S_estw) Xhx(1:length(P) - 1)];
    Cw = Cw + mu * err(k) * Xhx;       
end

% the result
figure;
subplot(2, 1, 1);
plot([1:T], err);
ylabel('Amplitude');
xlabel('time k');
legend('Noise redundent');
subplot(2, 1, 2);
plot([1:T], Yd);
hold on;
plot([1:T], Yd-err, 'r:');
ylabel('Amplitude');
xlabel('time k');
legend('Noise signal', 'Control signal');
