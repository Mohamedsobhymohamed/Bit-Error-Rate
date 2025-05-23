clc; clear;

%% Parameters
numBits = 1e6;                          % Total number of bits to transmit
SNR_dB_range = 0:3:60;                  % SNR values in dB
SNR_linear = 10.^(SNR_dB_range/10);    % Convert SNR from dB to linear scale

% Initialize BER arrays
ber_OOK_manual = [];
ber_PRK_manual = [];
ber_BFSK_manual = [];

%% Generate random binary data
bitStream = randi([0 1], 1, numBits);

%% Manual Modulation

% On-Off Keying (OOK): 0 -> 0, 1 -> 1
signal_OOK = bitStream;

% Phase Reversal Keying (PRK): 0 -> -1, 1 -> +1
signal_PRK = 2 * bitStream - 1;

% Binary Frequency Shift Keying (BFSK): 0 -> 1, 1 -> 1i (orthogonal)
signal_BFSK = ones(1, numBits);
signal_BFSK(bitStream == 1) = 1i;

%% Loop over all SNR values 
snr_idx = 1;
while snr_idx <= length(SNR_linear)
    snr_val = SNR_linear(snr_idx);

    % --- OOK ---
    power_OOK = mean(signal_OOK.^2);
    noisePower = power_OOK / snr_val;
    noise = sqrt(noisePower/2) * randn(1, numBits);
    received_OOK = signal_OOK + noise;
    decided_OOK = received_OOK > 0.5;
    ber = sum(decided_OOK ~= bitStream) / numBits;
    ber_OOK_manual = [ber_OOK_manual ber];

    % --- PRK ---
    power_PRK = mean(signal_PRK.^2);
    noisePower = power_PRK / snr_val;
    noise = sqrt(noisePower/2) * randn(1, numBits);
    received_PRK = signal_PRK + noise;
    decided_PRK = received_PRK > 0;
    ber = sum(decided_PRK ~= bitStream) / numBits;
    ber_PRK_manual = [ber_PRK_manual ber];

    % --- BFSK ---
    power_BFSK = mean(abs(signal_BFSK).^2);
    noisePower = power_BFSK / snr_val;
    noise = sqrt(noisePower/2) * (randn(1, numBits) + 1i * randn(1, numBits));
    received_BFSK = signal_BFSK + noise;
    decided_BFSK = imag(received_BFSK) > real(received_BFSK);
    ber = sum(decided_BFSK ~= bitStream) / numBits;
    ber_BFSK_manual = [ber_BFSK_manual ber];

    % Move to next SNR value
    snr_idx = snr_idx + 1;
end

%% Plot BER for manual modulations 
figure;
semilogy(SNR_dB_range, ber_OOK_manual, '-', 'Color', [0.85 0 0], 'LineWidth', 2.5); hold on;
semilogy(SNR_dB_range, ber_PRK_manual, '-', 'Color', [0 0.3 1], 'LineWidth', 2.5);
semilogy(SNR_dB_range, ber_BFSK_manual, '-', 'Color', [0 0.6 0], 'LineWidth', 2.5);
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER of OOK, PRK, BFSK (Manual)');
legend('OOK', 'PRK', 'BFSK');
grid on;

%% Built-in MATLAB Verification
ber_OOK_builtin = [];
ber_PRK_builtin = [];
ber_BFSK_builtin = [];
M_bfsk = 2;

snr_idx = 1;
while snr_idx <= length(SNR_dB_range)
    snr_dB = SNR_dB_range(snr_idx);

    % --- OOK ---
    mod_OOK = pammod(bitStream, 2);  % Maps 0 → 0, 1 → 1
    rx_OOK = awgn(mod_OOK, snr_dB, 'measured');
    demod_OOK = pamdemod(rx_OOK, 2);
    ber = sum(demod_OOK ~= bitStream) / numBits;
    ber_OOK_builtin = [ber_OOK_builtin ber];

    % --- PRK ---
    mod_PRK = pammod(bitStream, 2, 0, 'gray');
    mod_PRK = 2 * mod_PRK - 1;  % Scale to [-1, +1]
    rx_PRK = awgn(mod_PRK, snr_dB, 'measured');
    demod_PRK = double(rx_PRK > 0);
    ber = sum(demod_PRK ~= bitStream) / numBits;
    ber_PRK_builtin = [ber_PRK_builtin ber];

    % --- BFSK ---
    mod_BFSK = fskmod(bitStream, M_bfsk, 1, 100, 1e3);
    rx_BFSK = awgn(mod_BFSK, snr_dB, 'measured');
    demod_BFSK = fskdemod(rx_BFSK, M_bfsk, 1, 100, 1e3);
    ber = sum(demod_BFSK ~= bitStream) / numBits;
    ber_BFSK_builtin = [ber_BFSK_builtin ber];

    % Move to the next SNR value
    snr_idx = snr_idx + 1;
end

% Add to plot
semilogy(SNR_dB_range, ber_OOK_builtin, 'r--', 'LineWidth', 1.2);
semilogy(SNR_dB_range, ber_PRK_builtin, 'b--', 'LineWidth', 1.2);
semilogy(SNR_dB_range, ber_BFSK_builtin, 'g--', 'LineWidth', 1.2);
legend('OOK (manual)', 'PRK (manual)', 'BFSK (manual)', ...
       'OOK (built-in)', 'PRK (built-in)', 'BFSK (built-in)');

%% Part II: M-ASK Simulation for M = 2, 4, 8
M_values = [2 4 8];
ber_MASK_all = cell(1, length(M_values));
spacing = 2;

for idx = 1:length(M_values)
    M = M_values(idx);
    bitsPerSymbol = log2(M);
    numBits_adj = floor(numBits / bitsPerSymbol) * bitsPerSymbol;
    bits = randi([0 1], 1, numBits_adj);

    % Reshape into symbols
    bitMatrix = reshape(bits, bitsPerSymbol, []).';
    symbolIndices = bi2de(bitMatrix, 'left-msb');

    % Define symbol levels
    levels = ((0:M-1) - (M-1)/2) * spacing;
    symbols = levels(symbolIndices + 1);

    % BER vector for this M
    ber_vec = [];

    for snr_val = SNR_linear
        P = mean(symbols.^2);
        N0 = P / snr_val;
        noise = sqrt(N0/2) * randn(1, length(symbols));
        rx = symbols + noise;

        % Detection
        estimatedSymbols = zeros(1, length(rx));
        for k = 1:length(rx)
            [~, estimatedSymbols(k)] = min(abs(rx(k) - levels));
        end

        % Decode bits
        decodedBits = de2bi(estimatedSymbols - 1, bitsPerSymbol, 'left-msb');
        decodedBits = reshape(decodedBits.', 1, []);
        ber = sum(decodedBits ~= bits) / numBits_adj;
        ber_vec = [ber_vec ber];
    end

    ber_MASK_all{idx} = ber_vec;
end

%% Plot M-ASK BER curves 
figure;
semilogy(SNR_dB_range, ber_MASK_all{1}, '-o', 'Color', [0 0 0.8], 'LineWidth', 2.5); hold on;
semilogy(SNR_dB_range, ber_MASK_all{2}, '-s', 'Color', [0.9 0 0], 'LineWidth', 2.5);
semilogy(SNR_dB_range, ber_MASK_all{3}, '-^', 'Color', [0 0.6 0], 'LineWidth', 2.5);

xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER of M-ASK (M=2,4,8)');
legend('2-ASK', '4-ASK', '8-ASK');
grid on;

%% Part III: BPSK, QPSK, 4QAM Simulation
spacing = 2;
ber_BPSK = [];
ber_QPSK = [];
ber_QAM4 = [];

% --- BPSK ---
bits = randi([0 1], 1, numBits);
tx_BPSK = (2 * bits - 1) * (spacing/2);

for snr_val = SNR_linear
    P = mean(tx_BPSK.^2);
    N0 = P / snr_val;
    noise = sqrt(N0/2) * randn(1, length(tx_BPSK));
    rx = tx_BPSK + noise;
    decisions = double(rx > 0);
    ber = sum(decisions ~= bits) / numBits;
    ber_BPSK = [ber_BPSK ber];
end

% --- QPSK ---
bits = randi([0 1], 1, numBits);
I_bits = 2 * bits(1:2:end) - 1;
Q_bits = 2 * bits(2:2:end) - 1;
tx_QPSK = (I_bits + 1i * Q_bits) * (spacing/2) / sqrt(2);

for snr_val = SNR_linear
    P = mean(abs(tx_QPSK).^2);
    N0 = P / snr_val;
    noise = sqrt(N0/2) * (randn(1, length(tx_QPSK)) + 1i * randn(1, length(tx_QPSK)));
    rx = tx_QPSK + noise;
    dec_I = double(real(rx) > 0);
    dec_Q = double(imag(rx) > 0);
    estimatedBits = zeros(1, numBits);
    estimatedBits(1:2:end) = dec_I;
    estimatedBits(2:2:end) = dec_Q;
    ber = sum(estimatedBits ~= bits) / numBits;
    ber_QPSK = [ber_QPSK ber];
end

% --- 4QAM (Gray-mapped) ---
bits = randi([0 1], 1, numBits);
bitPairs = reshape(bits, 2, []).';
symbolIndices = bi2de(bitPairs, 'left-msb');
constellation = [ -1-1i, -1+1i, 1+1i, 1-1i ] * (spacing/2) / sqrt(2);
tx_QAM4 = constellation(symbolIndices + 1);

for snr_val = SNR_linear
    P = mean(abs(tx_QAM4).^2);
    N0 = P / snr_val;
    noise = sqrt(N0/2) * (randn(1, length(tx_QAM4)) + 1i * randn(1, length(tx_QAM4)));
    rx = tx_QAM4 + noise;

    % Detection
    estimatedSymbols = zeros(1, length(rx));
    for k = 1:length(rx)
        [~, estimatedSymbols(k)] = min(abs(rx(k) - constellation));
    end

    % Convert to bits
    decodedBits = de2bi(estimatedSymbols - 1, 2, 'left-msb');
    decodedBits = reshape(decodedBits.', 1, []);
    ber = sum(decodedBits ~= bits) / numBits;
    ber_QAM4 = [ber_QAM4 ber];
end

%% Final Plot: BPSK vs QPSK vs 4QAM vs 4ASK
figure;
semilogy(SNR_dB_range, ber_BPSK, 'b-o', 'LineWidth', 1.5); hold on;
semilogy(SNR_dB_range, ber_QPSK, 'r-s', 'LineWidth', 1.5);
semilogy(SNR_dB_range, ber_QAM4, 'g-^', 'LineWidth', 1.5);
semilogy(SNR_dB_range, ber_MASK_all{2}, 'k--d', 'LineWidth', 1.5); % 4ASK

xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER: BPSK vs QPSK vs 4QAM vs 4ASK');
legend('BPSK', 'QPSK', '4QAM', '4ASK');
grid on;
