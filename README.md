# CSP to Parisi

This code calculates the average satisfying fraction of a random sparse CSP, using the Parisi formula.

## Usage

Get Fourier levels $[|f_{=0}|^2,...,|f_{=k}|^2]$ of a Boolean predicate $f$:

```
python -i get_fourier_levels.py
```

```
>>> k = 3                       # number of variables
>>> inp = '01111111'            # specifies 1=T, 0=F for a predicate
>>> get_truth_table(k, inp)     # to see the truth table
>>> get_fourier_levels(k, inp)  # outputs Fourier levels
```
