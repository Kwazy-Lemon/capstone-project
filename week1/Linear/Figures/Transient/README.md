# Transient Analysis Results

This folder contains the results from the Week 1 transient analysis prototype.

## Circuit

The circuit consists of:

- A 10 V step voltage source
- R1 = 1 kΩ
- R2 = 2 kΩ
- C = 100 nF

The input voltage changes from 0 V to 10 V at \(t = 0\).

## Transient MNA

The transient analysis is performed using the Modified Nodal Analysis (MNA) formulation with the **Backward Euler** method.

The capacitor is discretized using:

\[
i_C =
\frac{C}{\Delta t}
\left(v_C^n-v_C^{n-1}\right)
\]

This introduces a time-dependent conductance:

\[
G_C = \frac{C}{\Delta t}
\]

and requires the voltage from the previous time step when constructing the right-hand-side vector.

## Results

The MATLAB prototype produces the following final values:

- Final \(V_{in}\) = 10.0000 V
- Final \(V_{out}\) = 6.6667 V
- Final \(I_{Vs}\) = -3.3333 mA

The output voltage starts near 0 V and gradually approaches the steady-state value of approximately 6.667 V.
<img width="461" height="456" alt="image" src="https://github.com/user-attachments/assets/870a731b-1c56-4a44-af13-4821ce921626" />


### Transient Response

<img width="444" height="245" alt="image" src="https://github.com/user-attachments/assets/34ac9aab-744c-409c-8df6-8ad97da18e1b" />

### Input and Output Voltage

<img width="438" height="243" alt="image" src="https://github.com/user-attachments/assets/87c06c28-4734-4b81-8aa9-930afc1caabc" />


## Summary

The results demonstrate the basic transient behavior of the RC circuit and verify the use of time stepping and previous time-step information in the transient MNA formulation.
