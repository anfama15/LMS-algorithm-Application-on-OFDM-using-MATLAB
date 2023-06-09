figure(2);
stem(abs(hf));
xlabel('Subcarrier label');
title('Subchannel gain');
% Plot the subchannel constellation scattering after OFDM
figure(3);
subplot(221);plot(p_data(:,1:800),'.') % subchannel 1 output
ylabel('Imaginary');
title('(a) Mixed OFDM without I/Q imbalance ');axis('square');
subplot(222);plot(n_data(:,1:800),'.'); % subchannel 10 output
ylabel('Imaginary');
title('(b) Mixed OFDM with I/Q imbalance ');axis('square');
subplot(223);plot(z_data(1,1:800),'.'); % subchannel 15 output
xlabel('Real');ylabel('Imaginary');
title('(c) Filtered Subchannel 1 output');axis('square');
subplot(224);plot(z_data(:,1:800),'b.'); % mixed subchannel output
xlabel('Real');ylabel('Imaginary');
title('(d) Mixed OFDM output');axis('square');
% Plot the average OFDM SER versus SER under "ideal channel"
% By Disabling 5 poor subcarriers, average SER can be reduced.
figure(4);
figc=semilogy(SNRa,Q,'k-',SNR,mean(SEReq),'b-o',...
SNR,mean([SEReq(1:14,:);SEReq(20:32,:)]),'b-s');
set(figc,'LineWidth',2);
legend('Ideal channel','Using all subcarriers','Disabling 5 poor subcarriers')
title('Average OFDM SER');
axis([1 30 1.e-4 1]);hold off;
xlabel('SNR (dB)');ylabel('Symbol Error Rate (SER)');
figure(5);
constellation(p_data)
