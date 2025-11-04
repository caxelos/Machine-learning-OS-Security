# PART A: FILTER DESIGN

We wish to design a discrete-time lowpass filter with the following specifications:

- **Passband** defined by frequency ωₚ = 0.10π  
- **Stopband** defined by frequency ωₛ = 0.30π  
- **Maximum ripple in the passband** defined by:

  \[
  R_p = -20 \log_{10} \frac{1 - \delta_1}{1 + \delta_1} = 1.00 \text{ dB}
  \]

- **Maximum ripple in the stopband** defined by:

  \[
  A_s = -20 \log_{10} \frac{\delta_2}{1 + \delta_1} = 40.00 \text{ dB}
  \]

---

## Tasks

### A1. 
Design the desired filter as an **FIR filter** using a **Hamming window**.

### A2. 
Also design the filter as an **IIR Butterworth filter** using the **bilinear transform**.

---

## For Each Case, Design and Plot:

- The **impulse response** of the filter  
- The **step response** of the filter  
- The **frequency response magnitude** in a logarithmic (dB) scale  
- The **group delay**  
- The **zero-pole diagram** (only for the IIR filter)

---

## Hint

Useful MATLAB commands include (among others):

```
fir1, sinc, hamming, buttord, butter, freqz, impz, stepz, grpdelay, zplane, tf, ztf, stem, figure, print
```
