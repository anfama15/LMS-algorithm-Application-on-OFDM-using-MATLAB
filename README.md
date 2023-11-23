# LMS-algorithm-Application-on-OFDM-using-MATLAB
The code is an implementation of an OFDM system with a polar signaling scheme and the application of an LMS filter for equalization.

The code starts with clearing the workspace and the figures using the clear and clf commands. Then, it defines the total number of data symbols L and the number of data frames Lfr. Next, it generates random signal data for polar signaling using 16QAM modulation. The generated signal data is stored in the variable s_data.

After generating the signal data, the code defines a channel in the time domain using the channel variable. The frequency response of the channel is obtained by computing its FFT using the fft function with a window size of 32. The data is then serialized and parallelized using the reshape function and the first data is equalized using an inverse of the frequency response of the channel. The data is then converted to the time domain using the ifft function.

Next, a cyclic prefix is added to the data to prevent inter-symbol interference. The power of the input signal to the channel is then set to 10/32 using the Psig variable. The channel output signal is obtained by filtering the input signal with the channel using the filter function.

The code then adds noise to the channel output signal to simulate channel noise. The alpha and beta variables are calculated to simulate I/Q imbalance in the system. The signal with I/Q imbalance is stored in the y_bar variable.

After adding the noise and simulating I/Q imbalance, the code equalizes the signal using the LMS filter. The mu variable is set to a convergence factor of 0.00018, and the filter is designed to have 32 taps. The LMS filter is implemented using a for loop, and the error signal en is computed at each iteration. The weights of the filter are updated using the iterative formula for the LMS algorithm. The equalized data is then converted to the frequency domain using the fft function and the frequency response of the channel is equalized.

Finally, the decision on the 16QAM data is made by comparing the equalized data with the original data. The symbol error rate (SER) is calculated by dividing the number of symbol errors by the total number of symbols. The SER is computed for different signal-to-noise ratios (SNRs) ranging from 0 to 30 dB. The code then computes the theoretical bit error rate (BER) using the erfc function and plots the BER and SER as a function of SNR. The ofdmaz function is then called to display the constellation plot of the received and equalized signals.

In summary, the code generates a random signal data for polar signaling and simulates a channel with noise and I/Q imbalance. The LMS filter is then used to equalize the signal and the symbol error rate is calculated for different SNRs. Finally, the BER and SER are plotted as a function of SNR, and the constellation plot of the received and equalized signals is displayed.


ofdmaz.m
This code generates several plots to visualize the performance of an OFDM system with an LMS filter applied. Here is a breakdown of what each part of the code does:

