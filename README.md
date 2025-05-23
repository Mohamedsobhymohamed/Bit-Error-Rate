This MATLAB project simulates and analyzes the **Bit Error Rate (BER)** performance of various digital modulation schemes under Additive White Gaussian Noise (AWGN) conditions. It covers manual and built-in implementations of modulation techniques and visualizes their performance using BER vs. SNR plots.

---

## 🧪 Modulation Schemes Simulated

### Part I: Manual vs Built-in BER Simulation
- **On-Off Keying (OOK)**
- **Phase Reversal Keying (PRK)** (equivalent to BPSK)
- **Binary Frequency Shift Keying (BFSK)**

### Part II: M-ASK Simulation
- M = 2, 4, 8

### Part III: Complex Modulations
- **BPSK**
- **QPSK**
- **4-QAM (Gray-mapped)**
- **Comparison with 4-ASK**

---

## 📁 File Overview

- `DigitalFinalProject.m`: Main script that runs the simulations, generates BER values, and produces comparative plots.

---

## 🚀 How to Run

1. Open MATLAB.
2. Load `DigitalFinalProject.m` into the editor.
3. Run the script using the `Run` button or `F5`.

> Ensure the Communications Toolbox is installed for built-in modulation functions (`pammod`, `fskmod`, `awgn`, etc.).

---

## 📊 Output

The script produces several **semilog BER vs. SNR plots** including:
- Manual vs built-in BER for OOK, PRK, BFSK
- BER for 2-ASK, 4-ASK, 8-ASK
- BER comparison of BPSK, QPSK, 4-QAM, and 4-ASK

---

## ⚙️ Parameters

| Parameter         | Value/Range              | Description                        |
|------------------|--------------------------|------------------------------------|
| `numBits`         | 1,000,000 (default)       | Number of bits to simulate         |
| `SNR_dB_range`    | 0:3:60                    | Range of SNR values in dB          |
| `spacing`         | 2                         | Symbol spacing for M-ASK/QAM       |

---

## 📚 Concepts Covered

- Digital modulation theory and practice
- BER computation under AWGN
- Modulation order impact (M-ary systems)
- Comparison between modulation techniques
- Gray coding and its BER benefits

---

## 📝 License

This project is intended for educational purposes. No commercial license is granted.

---

## 🙌 Acknowledgments

- MATLAB Communications Toolbox
- Standard definitions and formulas from digital communications theory
