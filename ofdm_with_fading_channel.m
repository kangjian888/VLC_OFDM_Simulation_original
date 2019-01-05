%OFDM Simulation in MATLAB
%The intial code is from Jordan Street's work. 
close all;
clc
clear all;
%% Simulation Parameters
%Moduluation method: BPSK, QPSK, 8PSK, 16QAM, 32QAM, 64QAM,
%Assume we use 14 channels of 256. 1,2 and 3 we use 128QAM, 4 and 5 we use 64QAM, 6,7,8 we use 32QAM
mod_methods = {'BPSK', 'QPSK', '8QAM', '16QAM', '32QAM', '64QAM','128QAM','256QAM'};
%Channel Information we use
Sub_channel_num_used = 80;
Channel_bit = [10,10,10,10,10,10,10,10];%Channel number of each modulation scheme
Total_bit_num_in_one_ofdm_symbol = 8*Channel_bit(8) +  7*Channel_bit(7) +  6*Channel_bit(6) +  5*Channel_bit(5) +  4*Channel_bit(4) +  3*Channel_bit(3) +  2*Channel_bit(2) +  1*Channel_bit(1);
%fft size
nfft = 512;
%size of cycle prefix extension
n_cpe = 16;

%Channel estimation method
ch_est_method = 'LS';
%ch_est_method = 'none';

%SNR in dB
snr = 10;%The unit is Db, the ratio of average power of signal and noise
%Number of channel taps
n_taps  = 8;%this is the same as the example, I think it is the pulse response
%Check whether the parameter is right or not?
if Sub_channel_num_used ~= sum(Channel_bit)
	disp('error: The summary of bit allocated channel is not equal to total used channel');
	exit;
end
if Sub_channel_num_used >= nfft/2 - 1%dc cannot be used in the simulation
	disp('error: The number of used subchannel is larger than fft number');
	exit;
end
%%Read the image and convert to binary format, we use one picture as our signal source
im =  imread('baboon.png');
im_bin_1 = im(:);
im_bin_2 = dec2bin(im_bin_1,8)';
im_bin = im_bin_2(:); %transfer the matrix to bit sequence
%% Transmisstor
%Binary Stream to Symbol
%We assume we just use 20 subchannel, the first five we use first five channel as 
%Step 1: Padding source signal
sys_rem = mod(Total_bit_num_in_one_ofdm_symbol-mod(length(im_bin),Total_bit_num_in_one_ofdm_symbol),Total_bit_num_in_one_ofdm_symbol);%calculate the bit not enought for one OFDM frame
padding = repmat('0',sys_rem,1);
im_bin_padded = [im_bin; padding];%combining in row
ofdm_symbol_total_num = length(im_bin_padded)/Total_bit_num_in_one_ofdm_symbol;
cons_data = reshape(im_bin_padded,Total_bit_num_in_one_ofdm_symbol,ofdm_symbol_total_num);
%Mapping source bit to symbol
%Step 2: Devide the signal source into different groups according modulation scheme, offline processing method
cons_data_1bit = cons_data(1:1*Channel_bit(1),:);
cons_data_2bit = cons_data(1*Channel_bit(1)+1:1*Channel_bit(1)+2*Channel_bit(2),:);
cons_data_3bit = cons_data(1*Channel_bit(1)+2*Channel_bit(2)+1:1*Channel_bit(1)+2*Channel_bit(2)+3*Channel_bit(3),:);
cons_data_4bit = cons_data(1*Channel_bit(1)+2*Channel_bit(2)+3*Channel_bit(3)+1:1*Channel_bit(1)+2*Channel_bit(2)+3*Channel_bit(3)+4*Channel_bit(4),:);
cons_data_5bit = cons_data(1*Channel_bit(1)+2*Channel_bit(2)+3*Channel_bit(3)+4*Channel_bit(4)+1:1*Channel_bit(1)+2*Channel_bit(2)+3*Channel_bit(3)+4*Channel_bit(4)+5*Channel_bit(5),:);
cons_data_6bit = cons_data(1*Channel_bit(1)+2*Channel_bit(2)+3*Channel_bit(3)+4*Channel_bit(4)+5*Channel_bit(5)+1:1*Channel_bit(1)+2*Channel_bit(2)+3*Channel_bit(3)+4*Channel_bit(4)+5*Channel_bit(5)+6*Channel_bit(6),:);
cons_data_7bit = cons_data(1*Channel_bit(1)+2*Channel_bit(2)+3*Channel_bit(3)+4*Channel_bit(4)+5*Channel_bit(5)+6*Channel_bit(6)+1:1*Channel_bit(1)+2*Channel_bit(2)+3*Channel_bit(3)+4*Channel_bit(4)+5*Channel_bit(5)+6*Channel_bit(6)+7*Channel_bit(7),:);
cons_data_8bit = cons_data(1*Channel_bit(1)+2*Channel_bit(2)+3*Channel_bit(3)+4*Channel_bit(4)+5*Channel_bit(5)+6*Channel_bit(6)+7*Channel_bit(7)+1:1*Channel_bit(1)+2*Channel_bit(2)+3*Channel_bit(3)+4*Channel_bit(4)+5*Channel_bit(5)+6*Channel_bit(6)+7*Channel_bit(7)+8*Channel_bit(8),:);
%Step 3: reshaping the each group for symbol mapping, and do the tranposition operation for mapping
cons_data_1bit_reshape = reshape(cons_data_1bit,1,ofdm_symbol_total_num*Channel_bit(1))';
cons_data_2bit_reshape = reshape(cons_data_2bit,2,ofdm_symbol_total_num*Channel_bit(2))';
cons_data_3bit_reshape = reshape(cons_data_3bit,3,ofdm_symbol_total_num*Channel_bit(3))';
cons_data_4bit_reshape = reshape(cons_data_4bit,4,ofdm_symbol_total_num*Channel_bit(4))';
cons_data_5bit_reshape = reshape(cons_data_5bit,5,ofdm_symbol_total_num*Channel_bit(5))';
cons_data_6bit_reshape = reshape(cons_data_6bit,6,ofdm_symbol_total_num*Channel_bit(6))';
cons_data_7bit_reshape = reshape(cons_data_7bit,7,ofdm_symbol_total_num*Channel_bit(7))';
cons_data_8bit_reshape = reshape(cons_data_8bit,8,ofdm_symbol_total_num*Channel_bit(8))';
%Step 4: transform from bit to dec for symbol mapping
cons_data_1bit_sys_id = bin2dec(cons_data_1bit_reshape);
cons_data_2bit_sys_id = bin2dec(cons_data_2bit_reshape);
cons_data_3bit_sys_id = bin2dec(cons_data_3bit_reshape);
cons_data_4bit_sys_id = bin2dec(cons_data_4bit_reshape);
cons_data_5bit_sys_id = bin2dec(cons_data_5bit_reshape);
cons_data_6bit_sys_id = bin2dec(cons_data_6bit_reshape);
cons_data_7bit_sys_id = bin2dec(cons_data_7bit_reshape);
cons_data_8bit_sys_id = bin2dec(cons_data_8bit_reshape);


%Adding this step to calculate the symbol error rate in each subchannel
cons_data_1bit_sys_id_reshape = reshape(cons_data_1bit_sys_id,Channel_bit(1),ofdm_symbol_total_num);
cons_data_2bit_sys_id_reshape = reshape(cons_data_2bit_sys_id,Channel_bit(2),ofdm_symbol_total_num);
cons_data_3bit_sys_id_reshape = reshape(cons_data_3bit_sys_id,Channel_bit(3),ofdm_symbol_total_num);
cons_data_4bit_sys_id_reshape = reshape(cons_data_4bit_sys_id,Channel_bit(4),ofdm_symbol_total_num);
cons_data_5bit_sys_id_reshape = reshape(cons_data_5bit_sys_id,Channel_bit(5),ofdm_symbol_total_num);
cons_data_6bit_sys_id_reshape = reshape(cons_data_6bit_sys_id,Channel_bit(6),ofdm_symbol_total_num);
cons_data_7bit_sys_id_reshape = reshape(cons_data_7bit_sys_id,Channel_bit(7),ofdm_symbol_total_num);
cons_data_8bit_sys_id_reshape = reshape(cons_data_8bit_sys_id,Channel_bit(8),ofdm_symbol_total_num);

%Step 4: Mapping and judging whether the matrix is empty or not
not_empty_flag = zeros(8,1);
if Channel_bit(1)~=0
	cons_data_1bit_symbol = qammod(cons_data_1bit_sys_id,2^1,'UnitAveragePower',true);
	cons_data_1bit_symbol_reshape = reshape(cons_data_1bit_symbol,Channel_bit(1),ofdm_symbol_total_num);
	not_empty_flag(1) = 1;
else
	cons_data_1bit_symbol_reshape = [];
end
if Channel_bit(2)~=0
	cons_data_2bit_symbol = qammod(cons_data_2bit_sys_id,2^2,'UnitAveragePower',true);
	cons_data_2bit_symbol_reshape = reshape(cons_data_2bit_symbol,Channel_bit(2),ofdm_symbol_total_num);
	not_empty_flag(2) = 1;
else
	cons_data_2bit_symbol_reshape = [];
end
if Channel_bit(3)~=0
	cons_data_3bit_symbol = qammod(cons_data_3bit_sys_id,2^3,'UnitAveragePower',true);
	cons_data_3bit_symbol_reshape = reshape(cons_data_3bit_symbol,Channel_bit(3),ofdm_symbol_total_num);
	not_empty_flag(3) = 1;
else
	cons_data_3bit_symbol_reshape = [];
end
if Channel_bit(4)~=0
	cons_data_4bit_symbol = qammod(cons_data_4bit_sys_id,2^4,'UnitAveragePower',true);
	cons_data_4bit_symbol_reshape = reshape(cons_data_4bit_symbol,Channel_bit(4),ofdm_symbol_total_num);
	not_empty_flag(4) = 1;
else
	cons_data_4bit_symbol_reshape = [];
end
if Channel_bit(5)~=0
	cons_data_5bit_symbol = qammod(cons_data_5bit_sys_id,2^5,'UnitAveragePower',true);
	cons_data_5bit_symbol_reshape = reshape(cons_data_5bit_symbol,Channel_bit(5),ofdm_symbol_total_num);
	not_empty_flag(5) = 1;
else
	cons_data_5bit_symbol_reshape = [];
end
if Channel_bit(6)~=0
	cons_data_6bit_symbol = qammod(cons_data_6bit_sys_id,2^6,'UnitAveragePower',true);
	cons_data_6bit_symbol_reshape = reshape(cons_data_6bit_symbol,Channel_bit(6),ofdm_symbol_total_num);
	not_empty_flag(6) = 1;
else
	cons_data_6bit_symbol_reshape = [];
end
if Channel_bit(7)~=0
	cons_data_7bit_symbol = qammod(cons_data_7bit_sys_id,2^7,'UnitAveragePower',true);
	cons_data_7bit_symbol_reshape = reshape(cons_data_7bit_symbol,Channel_bit(7),ofdm_symbol_total_num);
	not_empty_flag(7) = 1;
else
	cons_data_7bit_symbol_reshape = [];
end
if Channel_bit(8)~=0
	cons_data_8bit_symbol = qammod(cons_data_8bit_sys_id,2^8,'UnitAveragePower',true);
	cons_data_8bit_symbol_reshape = reshape(cons_data_8bit_symbol,Channel_bit(8),ofdm_symbol_total_num);
	not_empty_flag(8) = 1;
else
	cons_data_8bit_symbol_reshape = [];
end
%Step 5: Combining all matrix together
cons_data_bit_symbol_reshape = [cons_data_8bit_symbol_reshape;cons_data_7bit_symbol_reshape;cons_data_6bit_symbol_reshape;cons_data_5bit_symbol_reshape;cons_data_4bit_symbol_reshape;cons_data_3bit_symbol_reshape;cons_data_2bit_symbol_reshape;cons_data_1bit_symbol_reshape];

%Step 6: Mapping singal in fft matrix, do not use dc and N+1 subchannel, to make the singal real, we use hermistian matrix.
%C_k = C*_(2N-k+2)
cons_data_bit_symbol_reshape_conj = conj(cons_data_bit_symbol_reshape);
carrier = 2:(Sub_channel_num_used+1);
conjugate_carrier = nfft - carrier + 2;
ofdm_modulation_matrix = zeros(nfft,ofdm_symbol_total_num);
ofdm_modulation_matrix(carrier,:) = cons_data_bit_symbol_reshape;
ofdm_modulation_matrix(conjugate_carrier,:) = cons_data_bit_symbol_reshape_conj;
%Step 7: Do IFFT Modulation
Sending_Sequence = ifft(ofdm_modulation_matrix,nfft,1);%to row to do the ifft transformation
%Step 8: Adding CP
Sending_Sequence_CP = [Sending_Sequence(end-n_cpe+1:end,:);Sending_Sequence];
Sending_Sequence_CP_reshape = Sending_Sequence_CP(:);
%% Cearteing a channel
%AWGN Part
data_pwr = mean(abs(Sending_Sequence_CP_reshape.^2));
noise_pwd = data_pwr/10^(snr/10); %we could calculate the noise power.
Sending_Sequence_added_whitenoise = awgn(Sending_Sequence_CP_reshape,snr,'measured');
%Apply fading channel
g = exp(-(0:n_taps-1));
g = g/norm(g);
Sending_Sequence_after_fading_channel = conv(Sending_Sequence_added_whitenoise,g,'same');%The size of output signal is the same with imput

%%Receiver Part Simulation
%Step 1: remove the prefix extension and shoft from serial to parallel
Recieving_Sequence = Sending_Sequence_after_fading_channel;
%Recieving_Sequence = Sending_Sequence_CP_reshape;%this is just for test the system, transmitter signal to receiver directly
Recieving_Sequence_reshape = reshape(Recieving_Sequence,nfft+n_cpe,length(Recieving_Sequence)/(nfft+n_cpe));
Recieving_Sequence_reshape_remove_cp = Recieving_Sequence_reshape(n_cpe+1:end,:);

%Step 2: Move to the frequency domain
Recieving_hat_Sequence = fft(Recieving_Sequence_reshape_remove_cp,nfft,1);

%Step 3: Channel Estimation
if n_taps > 1
	switch (ch_est_method)
		case 'none'
		case 'LS'
			G = Recieving_hat_Sequence(:,1)./ofdm_modulation_matrix(:,1);
			Recieving_hat_Sequence = Recieving_hat_Sequence./repmat(G,1,size(Recieving_hat_Sequence,2));
	end
end

%Step 4: Extracting useful information
rece_carrier = 2:(Sub_channel_num_used+1);
Recieving_hat_symbol = Recieving_hat_Sequence(rece_carrier,:);

%Step 5: Demodulation bit information from modulated signal
%Allocation table should be generated by receiver and transmit to transmitter through feedback path
%The demodulation process is the inverse process of modulation process
rece_not_empty_flag = zeros(8,1);
if Channel_bit(1)~=0
	rece_data_1bit_symbol = Recieving_hat_symbol(end-Channel_bit(1)+1:end,:);
	rece_data_1bit_symbol_reshape = reshape(rece_data_1bit_symbol,Channel_bit(1)*ofdm_symbol_total_num,1);
	rece_data_1bit_sys_id = qamdemod(rece_data_1bit_symbol_reshape,2^1,'UnitAveragePower',true);
	rece_data_1bit_sys_id_reshape = reshape(rece_data_1bit_sys_id,Channel_bit(1),ofdm_symbol_total_num);
	rece_data_1bit = dec2bin(rece_data_1bit_sys_id);
	rece_data_1bit_reshape = reshape(rece_data_1bit',1*Channel_bit(1),ofdm_symbol_total_num);
	rece_not_empty_flag(1) = 1;
else
	rece_data_1bit_reshape = [];
end
if Channel_bit(2)~=0
	rece_data_2bit_symbol = Recieving_hat_symbol(end-Channel_bit(1)-Channel_bit(2)+1:end-Channel_bit(1),:);
	rece_data_2bit_symbol_reshape = reshape(rece_data_2bit_symbol,Channel_bit(2)*ofdm_symbol_total_num,1);
	rece_data_2bit_sys_id = qamdemod(rece_data_2bit_symbol_reshape,2^2,'UnitAveragePower',true);
	rece_data_2bit_sys_id_reshape = reshape(rece_data_2bit_sys_id,Channel_bit(2),ofdm_symbol_total_num);
	rece_data_2bit = dec2bin(rece_data_2bit_sys_id);
	rece_data_2bit_reshape = reshape(rece_data_2bit',2*Channel_bit(2),ofdm_symbol_total_num);
	rece_not_empty_flag(2) = 1;
else
	rece_data_2bit_reshape= [];
end
if Channel_bit(3)~=0
	rece_data_3bit_symbol = Recieving_hat_symbol(end-Channel_bit(1)-Channel_bit(2)-Channel_bit(3)+1:end-Channel_bit(1)-Channel_bit(2),:);
	rece_data_3bit_symbol_reshape = reshape(rece_data_3bit_symbol,Channel_bit(3)*ofdm_symbol_total_num,1);
	rece_data_3bit_sys_id = qamdemod(rece_data_3bit_symbol_reshape,2^3,'UnitAveragePower',true);
	rece_data_3bit_sys_id_reshape = reshape(rece_data_3bit_sys_id,Channel_bit(3),ofdm_symbol_total_num);
	rece_data_3bit = dec2bin(rece_data_3bit_sys_id);
	rece_data_3bit_reshape = reshape(rece_data_3bit',3*Channel_bit(3),ofdm_symbol_total_num);
	rece_not_empty_flag(3) = 1;
else
	rece_data_3bit_reshape = [];
end
if Channel_bit(4)~=0
	rece_data_4bit_symbol = Recieving_hat_symbol(end-Channel_bit(1)-Channel_bit(2)-Channel_bit(3)-Channel_bit(4)+1:end-Channel_bit(1)-Channel_bit(2)-Channel_bit(3),:);
	rece_data_4bit_symbol_reshape = reshape(rece_data_4bit_symbol,Channel_bit(4)*ofdm_symbol_total_num,1);
	rece_data_4bit_sys_id = qamdemod(rece_data_4bit_symbol_reshape,2^4,'UnitAveragePower',true);
	rece_data_4bit_sys_id_reshape = reshape(rece_data_4bit_sys_id,Channel_bit(4),ofdm_symbol_total_num);
	rece_data_4bit = dec2bin(rece_data_4bit_sys_id);
	rece_data_4bit_reshape = reshape(rece_data_4bit',4*Channel_bit(4),ofdm_symbol_total_num);
	rece_not_empty_flag(4) = 1;
else
	rece_data_4bit_reshape = [];
end
if Channel_bit(5)~=0
	rece_data_5bit_symbol = Recieving_hat_symbol(end-Channel_bit(1)-Channel_bit(2)-Channel_bit(3)-Channel_bit(4)-Channel_bit(5)+1:end-Channel_bit(1)-Channel_bit(2)-Channel_bit(3)-Channel_bit(4),:);
	rece_data_5bit_symbol_reshape = reshape(rece_data_5bit_symbol,Channel_bit(5)*ofdm_symbol_total_num,1);
	rece_data_5bit_sys_id = qamdemod(rece_data_5bit_symbol_reshape,2^5,'UnitAveragePower',true);
	rece_data_5bit_sys_id_reshape = reshape(rece_data_5bit_sys_id,Channel_bit(5),ofdm_symbol_total_num);
	rece_data_5bit = dec2bin(rece_data_5bit_sys_id);
	rece_data_5bit_reshape = reshape(rece_data_5bit',5*Channel_bit(5),ofdm_symbol_total_num);
	rece_not_empty_flag(5) = 1;
else
	rece_data_5bit_reshape = [];
end
if Channel_bit(6)~=0
	rece_data_6bit_symbol = Recieving_hat_symbol(end-Channel_bit(1)-Channel_bit(2)-Channel_bit(3)-Channel_bit(4)-Channel_bit(5)-Channel_bit(6)+1:end-Channel_bit(1)-Channel_bit(2)-Channel_bit(3)-Channel_bit(4)-Channel_bit(5),:);
	rece_data_6bit_symbol_reshape = reshape(rece_data_6bit_symbol,Channel_bit(6)*ofdm_symbol_total_num,1);
	rece_data_6bit_sys_id = qamdemod(rece_data_6bit_symbol_reshape,2^6,'UnitAveragePower',true);
	rece_data_6bit_sys_id_reshape = reshape(rece_data_6bit_sys_id,Channel_bit(6),ofdm_symbol_total_num);
	rece_data_6bit = dec2bin(rece_data_6bit_sys_id);
	rece_data_6bit_reshape = reshape(rece_data_6bit',6*Channel_bit(6),ofdm_symbol_total_num);
	rece_not_empty_flag(6) = 1;
else
	rece_data_6bit_reshape = [];
end
if Channel_bit(7)~=0
	rece_data_7bit_symbol = Recieving_hat_symbol(end-Channel_bit(1)-Channel_bit(2)-Channel_bit(3)-Channel_bit(4)-Channel_bit(5)-Channel_bit(6)-Channel_bit(7)+1:end-Channel_bit(1)-Channel_bit(2)-Channel_bit(3)-Channel_bit(4)-Channel_bit(5)-Channel_bit(6),:);
	rece_data_7bit_symbol_reshape = reshape(rece_data_7bit_symbol,Channel_bit(7)*ofdm_symbol_total_num,1);
	rece_data_7bit_sys_id = qamdemod(rece_data_7bit_symbol_reshape,2^7,'UnitAveragePower',true);
	rece_data_7bit_sys_id_reshape = reshape(rece_data_7bit_sys_id,Channel_bit(7),ofdm_symbol_total_num);
	rece_data_7bit = dec2bin(rece_data_7bit_sys_id);
	rece_data_7bit_reshape = reshape(rece_data_7bit',7*Channel_bit(7),ofdm_symbol_total_num);
	rece_not_empty_flag(7) = 1;
else
	rece_data_7bit_reshape = [];
end
if Channel_bit(8)~=0
	rece_data_8bit_symbol = Recieving_hat_symbol(end-Channel_bit(1)-Channel_bit(2)-Channel_bit(3)-Channel_bit(4)-Channel_bit(5)-Channel_bit(6)-Channel_bit(7)-Channel_bit(8)+1:end-Channel_bit(1)-Channel_bit(2)-Channel_bit(3)-Channel_bit(4)-Channel_bit(5)-Channel_bit(6)-Channel_bit(7),:);
	rece_data_8bit_symbol_reshape = reshape(rece_data_8bit_symbol,Channel_bit(8)*ofdm_symbol_total_num,1);
	rece_data_8bit_sys_id = qamdemod(rece_data_8bit_symbol_reshape,2^8,'UnitAveragePower',true);
	rece_data_8bit_sys_id_reshape = reshape(rece_data_8bit_sys_id,Channel_bit(8),ofdm_symbol_total_num);
	rece_data_8bit = dec2bin(rece_data_8bit_sys_id);
	rece_data_8bit_reshape = reshape(rece_data_8bit',8*Channel_bit(8),ofdm_symbol_total_num);
	rece_not_empty_flag(8) = 1;
else
	rece_data_8bit_reshape = [];
end
rece_data_bit = [rece_data_1bit_reshape;rece_data_2bit_reshape;rece_data_3bit_reshape;rece_data_4bit_reshape;rece_data_5bit_reshape;rece_data_6bit_reshape;rece_data_7bit_reshape;rece_data_8bit_reshape];
%Step 6: Removing padding and recovery the image
rece_im_bin_1 = rece_data_bit(:);
rece_im_bin_2 = rece_im_bin_1(1:end-sys_rem);
rece_im_bin_3 = reshape(rece_im_bin_2,8,length(rece_im_bin_2)/8);
rece_im_bin_4 = uint8(bin2dec(rece_im_bin_3'));
rece_im_bin_5 = reshape(rece_im_bin_4,size(im));
%% Generate plot
% First original image
figure(1)
subplot(2,2,1);
imshow(im);
title('\bfTransmit Image');
subplot(2,2,2);
plot(Sending_Sequence_CP_reshape);
title('\bfSending singal sequence');
subplot(2,2,3);
plot(Sending_Sequence_CP(:,1));
title('\bfSending singal sequence one OFDM symbol');
subplot(2,2,4);
imshow(rece_im_bin_5);
title('\bfReceive Image');
% Sending Constellation image
figure(2)
title('\bfTransmission Singal Constellation')
if not_empty_flag(1)
	subplot(2,4,1);
	plot(cons_data_1bit_symbol,zeros(1,length(cons_data_1bit_symbol)),'x','linewidth',2,'markersize',10);
	title('\bfBPSK');
    xlim([-2 2]);
    ylim([-2 2]);
else
	subplot(2,4,1);
	xlim([-2 2]);
	ylim([-2 2]);
	title(sprintf('No data used in BPSK modulation scheme'));
end
if not_empty_flag(2)
	subplot(2,4,2);
	plot(cons_data_2bit_symbol,'x','linewidth',2,'markersize',10);
	title('\bfQPSK');
    xlim([-2 2]);
    ylim([-2 2]);
else
	subplot(2,4,2);
	xlim([-2 2]);
	ylim([-2 2]);	
	title(sprintf('No data used in QPSK modulation scheme'));
end
if not_empty_flag(3)
	subplot(2,4,3);
	plot(cons_data_3bit_symbol,'x','linewidth',2,'markersize',10);
	title('\bf8QAM');
    xlim([-2 2]);
    ylim([-2 2]);
else
	subplot(2,4,3);
	title(sprintf('No data used in 8QAM modulation scheme'));
	xlim([-2 2]);
	ylim([-2 2]);
end
if not_empty_flag(4)
	subplot(2,4,4);
	plot(cons_data_4bit_symbol,'x','linewidth',2,'markersize',10);
	title('\bf16QAM');
    xlim([-2 2]);
    ylim([-2 2]);
else
	subplot(2,4,4);
	xlim([-2 2]);
	ylim([-2 2]);
	title(sprintf('No data used in 16QAM modulation scheme'));
end
if not_empty_flag(5)
	subplot(2,4,5);
	plot(cons_data_5bit_symbol,'x','linewidth',2,'markersize',10);
	title('\bf32QAM');
    xlim([-2 2]);
    ylim([-2 2]);
else
	subplot(2,4,5);
	xlim([-2 2]);
	ylim([-2 2]);
	title(sprintf('No data used in 32QAM modulation scheme'));
end
if not_empty_flag(6)
	subplot(2,4,6);
	plot(cons_data_6bit_symbol,'x','linewidth',2,'markersize',10);
	title('\bf64QAM');
    xlim([-2 2]);
    ylim([-2 2]);
else
	subplot(2,4,6);
	xlim([-2 2]);
	ylim([-2 2]);
	title(sprintf('No data used in 64QAM modulation scheme'));
end
if not_empty_flag(7)
	subplot(2,4,7);
	plot(cons_data_7bit_symbol,'x','linewidth',2,'markersize',10);
	title('\bf128QAM');
    xlim([-2 2]);
    ylim([-2 2]);
else
	subplot(2,4,7);
	xlim([-2 2]);
	ylim([-2 2]);
	title(sprintf('No data used in 128QAM modulation scheme'));
end
if not_empty_flag(8)
	subplot(2,4,8);
	plot(cons_data_8bit_symbol,'x','linewidth',2,'markersize',10);
	title('\bf256QAM');
    xlim([-2 2]);
	ylim([-2 2]);
else
	subplot(2,4,8);
	xlim([-2 2]);
	ylim([-2 2]);
	title(sprintf('No data used in 256QAM modulation scheme'));
end

% Receiver Constellation image
figure(3)
title('\bfReceiver Singal Constellation')
if rece_not_empty_flag(1)
	subplot(2,4,1);
	plot(rece_data_1bit_symbol_reshape,'x','linewidth',2,'markersize',10);
	title('\bfBPSK');
    xlim([-2 2]);
    ylim([-2 2]);
else
	subplot(2,4,1);
	xlim([-2 2]);
	ylim([-2 2]);
	title(sprintf('No data used in BPSK modulation scheme'));
end
if rece_not_empty_flag(2)
	subplot(2,4,2);
	plot(rece_data_2bit_symbol_reshape,'x','linewidth',2,'markersize',10);
	title('\bfQPSK');
    xlim([-2 2]);
    ylim([-2 2]);
else
	subplot(2,4,2);
	xlim([-2 2]);
	ylim([-2 2]);	
	title(sprintf('No data used in QPSK modulation scheme'));
end
if rece_not_empty_flag(3)
	subplot(2,4,3);
	plot(rece_data_3bit_symbol_reshape,'x','linewidth',2,'markersize',10);
	title('\bf8QAM');
    xlim([-2 2]);
    ylim([-2 2]);
else
	subplot(2,4,3);
	title(sprintf('No data used in 8QAM modulation scheme'));
	xlim([-2 2]);
	ylim([-2 2]);
end
if rece_not_empty_flag(4)
	subplot(2,4,4);
	plot(rece_data_4bit_symbol_reshape,'x','linewidth',2,'markersize',10);
	title('\bf16QAM');
    xlim([-2 2]);
    ylim([-2 2]);
else
	subplot(2,4,4);
	xlim([-2 2]);
	ylim([-2 2]);
	title(sprintf('No data used in 16QAM modulation scheme'));
end
if rece_not_empty_flag(5)
	subplot(2,4,5);
	plot(rece_data_5bit_symbol_reshape,'x','linewidth',2,'markersize',10);
	title('\bf32QAM');
    xlim([-2 2]);
    ylim([-2 2]);
else
	subplot(2,4,5);
	xlim([-2 2]);
	ylim([-2 2]);
	title(sprintf('No data used in 32QAM modulation scheme'));
end
if rece_not_empty_flag(6)
	subplot(2,4,6);
	plot(rece_data_6bit_symbol_reshape,'x','linewidth',2,'markersize',10);
	title('\bf64QAM');
    xlim([-2 2]);
    ylim([-2 2]);
else
	subplot(2,4,6);
	xlim([-2 2]);
	ylim([-2 2]);
	title(sprintf('No data used in 64QAM modulation scheme'));
end
if rece_not_empty_flag(7)
	subplot(2,4,7);
	plot(rece_data_7bit_symbol_reshape,'x','linewidth',2,'markersize',10);
	title('\bf128QAM');
    xlim([-2 2]);
    ylim([-2 2]);
else
	subplot(2,4,7);
	xlim([-2 2]);
	ylim([-2 2]);
	title(sprintf('No data used in 128QAM modulation scheme'));
end
if rece_not_empty_flag(8)
	subplot(2,4,8);
	plot(rece_data_8bit_symbol_reshape,'x','linewidth',2,'markersize',10);
	title('\bf256QAM');
    xlim([-2 2]);
	ylim([-2 2]);
else
	subplot(2,4,8);
	xlim([-2 2]);
	ylim([-2 2]);
	title(sprintf('No data used in 256QAM modulation scheme'));
end

%%SER and BER Calculation
%Part 1: SER Per Channel Calculation
