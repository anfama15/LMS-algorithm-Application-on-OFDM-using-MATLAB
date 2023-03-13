% Clear the workspace and figure
clear; clf;

% Define parameters
L = 160000; % Total data symbols in experiment
Lfr = L/32; % Number of data frames
channel = [0.3, -0.5, 0, 1, 0.2, -0.3]; % Channel in time-domain

% Generate random signal data for polar signaling
s_data = 4*round(rand(L, 1)) + 2*round(rand(L, 1)) - 3 ...
        + 1i*(4*round(rand(L, 1)) + 2*round(rand(L, 1)) - 3); % 16QAM modulation

% Reshape signal data for S/P conversion
p_data = reshape(s_data, 32, Lfr);

% IDFT to convert to time-domain
p_td = ifft(p_data);

% Cyclic prefix
p_cyc = [p_td(end-4:end, :); p_td]; % Add cyclic prefix
s_cyc = reshape(p_cyc, 37*Lfr, 1); % P/S conversion

% Find the channel in frequency-domain
hf = fft(channel, 32);

% Generate channel output signal
chsout = filter(channel, 1, s_cyc);

% Add complex Gaussian noise
noiseq = (randn(37*Lfr, 1) + 1i*randn(37*Lfr, 1));

% Calculate average channel input power
Psig = 10/32;

% Define SNR range
SNR = 0:30;

% Define I/Q imbalance parameters
alpha = chsout * 0.3;
y_conj = conj(chsout);
beta = y_conj * 0.7;
y_bar = alpha + beta;

% Initialize SEReq array
SEReq = [];

% Apply LMS algorithm for each SNR level
for ii = 1:length(SNR)
    
    % Add noise to channel output signal
    Asig = sqrt(Psig * 10^(-SNR(ii)/10)) * norm(channel);
    y_t = chsout + Asig * noiseq;
    
    % Convert signal to parallel form for LMS equalization
    n_para = reshape(y_bar, 37, Lfr);
    n_disc = n_para(6:37, :);
    nhat_para = fft(n_disc);
    n_data = inv(diag(hf)) * nhat_para;
    
    % Calculate convergence factor for LMS
    mu = 0.00018; % mu convergence factor (step size)
    
    % Define LMS filter parameters
    M = 32; % Number of filter taps
    itr = length(n_data); % Number of iterations
    en = zeros(itr, 1); % Error sequence
    W = zeros(M, itr); % Weighting parameter matrix
    
    % Iterative calculation for LMS filter
    for k = M:itr
        x = n_data(k:-1:k-M+1);
        y = W(:, k-1).' * x;
        en(k) = p_data(k) - y;
        W(:, k) = W(:, k-1) + 2*mu*en(k)*x;
    end
    
    % Get output sequence of the filter with the last best estimate
    yn = inf * ones(size(n_data));
    for k = M:length(n_data)
        x = n_data(k:-1:k-M+1);
        yn(k) = W(:, end).'
end
%%
x_para=reshape(yn,37,Lfr); % S/P conversion
x_disc=x_para(6:37,:); % discard tails
xhat_para=fft(x_disc); % FFT back to f-domain
z_data=inv(diag(hf))*xhat_para; % f-domain equalizing
% compute the 16QAM decision after equalization
deq=sign(real(z_data))+sign(real(z_data)-2)+sign(real(z_data)+2)+...
1i*(sign(imag(z_data))+sign(imag(z_data)-2)+sign(imag(z_data)+2));
% Now compare against the original data to compute SER
SEReq=[SEReq sum(p_data~=deq,2)/Lfr];
end
for ii=1:9,
    SNRa(ii)=2*ii-2;
    Q(ii)=3*0.5*erfc(sqrt((2*10^(SNRa(ii)*0.1)/5)/2));
%Compute the Analytical BER
end
ofdmaz;
