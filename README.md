# CSP to Parisi

This code calculates the average satisfying fraction of a random sparse CSP, using the Parisi formula.

## Usage

Get Fourier weights $[|f_{=0}|^2,...,|f_{=k}|^2]$ of a Boolean predicate $f$:

```bash
python -i get_fourier_weights.py
```

```python
>>> k = 3                       # number of variables
>>> inp = '01111111'            # specifies 1=T, 0=F for a predicate
>>> get_truth_table(k, inp)     # to see the truth table
>>> get_fourier_weights(k, inp)  # outputs Fourier weights
```

#### Note
This assumes the predicate has codomain ${-1, +1}$. 
If you need ${0,1}$ output, divide all nonzero Fourier weights by 4. 
(The weight at zero is the square of the average value of the function.)
