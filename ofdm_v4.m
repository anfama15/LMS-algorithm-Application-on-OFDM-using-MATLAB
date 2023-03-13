
clear;clf;
L=160000; % Total data symbols in experiment is 
Lfr=L/32; % number of data frames
% Generating random signal data for polar signaling
s_data=4*round(rand(L,1))+2*round(rand( L,1))-3+...
    +1i*(4*round(rand(L,1))+2*round(rand(L,1))-3); %16QAM Modulation
abc=s_data;
channel=[0.3 -0.5 0 1 .2 -0.3]; % channel in t-domain
hf=fft(channel,32); % find the channel in f-domain
p_data=reshape(s_data,32,Lfr); % S/P conversion
first_data=inv(diag(hf))*p_data;
p_td=ifft(p_data); % IDFT to convert to t-domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%cyclic prefix%%%%%%%%%%%%%%%%%%%%%%%%%
p_cyc=[p_td(end-4:end,:);p_td]; % add cyclic prefix
s_cyc=reshape(p_cyc,37*Lfr,1); % P/S conversion
Psig=10/32; % average channel input power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chsout=filter(channel,1,s_cyc); % generate channel output signal
clear p_td p_cyc s_data s_cyc; % release some memory
noiseq=(randn(37*Lfr,1)+1i*randn(37*Lfr,1));
SEReq=[];
for ii=1:31,
SNR(ii)=ii-1; % SNR in dB
% Asig=sqrt(Psig*10^(-SNR(ii)/10))*norm(channel);
% y_t=chsout+Asig*noiseq; % Add noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%I/Q Imbalance%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=chsout*0.3;
y_conj=conj(chsout);
beta=y_conj*0.7;
y_bar=alpha+beta;
%%
%Noise
n_para=reshape(y_bar,37,Lfr); % S/P conversion
n_disc=n_para(6:37,:); % discard tails
nhat_para=fft(n_disc); % FFT back to f-domain
n_data=inv(diag(hf))*nhat_para; % f-domain equalizing

%%
%LMS
%%
%dereiving an algo to calculate convergence factor.
mat=corrcoef(y_bar);
xn=y_bar;
dn=chsout;
M=32;
auto=xcorr(xn);
%eigen=eig(auto);
%mu=max(eigen);
%Px=2/(32*auto);
%mu=max(Px);
%%

mu=0.00018; % mu convergence factor (step size) 
itr = length(xn); 
en = zeros(itr,1);     % error sequence, en(k) represents the error between the expected output and the actual input at the k-th iteration        
W  = zeros(M,itr);     % Each row represents a weighting parameter, each column represents - iterations, and the initial value is 0
% iterative calculation        
for k = M:itr          % k-th iteration        
    x = xn(k:-1:k-M+1); % input of filter M taps       
    y = W(:,k-1).' * x; % output of the filter        
    en(k) = dn(k) - y ; % error at iteration k       
    % Iterative formula for filter weight calculation
    W(:,k) = W(:,k-1) + 2*mu*en(k)*x;
end
% The output sequence of the filter when finding the optimal r If there is no yn return parameter, the following can be omitted
yn = inf * ones(size(xn)); % inf means infinity
for k = M:length(xn)
    x = xn(k:-1:k-M+1);
    yn(k) = W(:,end).'* x; % Get the output with the last best estimate
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
