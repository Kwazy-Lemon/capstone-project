# AC Analysis Results

This folder contains the results from the Week 1 AC analysis prototype.

## Circuit

The circuit consists of:

- A 10 V AC voltage source
- R1 = 1 kΩ
- R2 = 2 kΩ
- C = 100 nF

The circuit is analyzed over a frequency range from 1 Hz to 1 MHz.

## AC MNA

The AC analysis extends the MNA formulation into the frequency domain.

For the capacitor, the frequency-dependent admittance is:

\[
Y_C = j\omega C
\]

The resulting MNA system is solved at each frequency point to obtain the output voltage response.

## Results

The MATLAB prototype produces the following results:

- Low-frequency output magnitude: approximately 6.667 V
- Output magnitude decreases as frequency increases
- The output phase approaches -90° at high frequency
- At 1 MHz, the output magnitude is approximately 0.016 V
<img width="458" height="492" alt="image" src="https://github.com/user-attachments/assets/fbbbea6d-d639-45a3-a638-dd71af917190" />

### Magnitude Response

<img width="436" height="241" alt="image" src="https://github.com/user-attachments/assets/c325cb8b-cfff-48e8-a081-257e2cf589c0" />


### Phase Response

<img width="445" height="242" alt="image" src="https://github.com/user-attachments/assets/46ab2a6f-2f09-4ff5-86e5-d9dcd061fa95" />


## Summary

The results demonstrate the expected frequency-dependent behavior of the RC circuit and verify the basic AC MNA formulation used in the MATLAB prototype.
