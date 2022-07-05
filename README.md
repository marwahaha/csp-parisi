# CSP to Parisi

This Python code calculates the average satisfying fraction of a random sparse CSP, using the Parisi formula.

packages required: SymPy, NumPy, SciPy

## Usage

```bash
python -i csp_parisi.py
```

To get the maximum satisfying fraction of a Boolean predicate $f$:

```python
>>> k = 3             # number of variables
>>> inp = '01111111'  # specifies 1=T, 0=F for a predicate
>>> get_truth_table(k, inp)
>>> weights = get_fourier_weights(k, inp)
>>> value = parisi_one_jump_approx(weights)
```

How this works is defined below:


### Get Fourier weights $[|f_{=0}|^2,...,|f_{=k}|^2]$ of a Boolean predicate $f$


```python
>>> k = 3                       # number of variables
>>> inp = '01111111'            # specifies 1=T, 0=F for a predicate
>>> get_truth_table(k, inp)     # to see the truth table
>>> get_fourier_weights(k, inp) # outputs Fourier weights
```

**Note**:
This assumes the predicate has codomain ${-1, +1}$.
If you need ${0,1}$ output, divide all nonzero Fourier weights by 4.
(The weight at zero is the square of the average value of the function.)


### Get approximate Parisi value, given list of mixture function coefficients $[c_0^2, c_1^2,...]$

```bash
python -i get_parisi_value.py
```

```python
C_p_squared = [0, 0, 0, 1/2]
>>> parisi_one_jump_approx(C_p_squared) # outputs Parisi value
```

This code only optimizes over piecewise constant functions with one jump.
It is reasonably fast (~10 seconds to run), and is >99.5% accurate on known results.
You can run `verify_parisi_one_jump_accuracy()` to double-check the accuracy claims.

### Get more accurate Parisi value, given list of mixture function coefficients $[c_0^2, c_1^2,...]$

```bash
python -i get_parisi_value.py
```

```python
C_p_squared = [0, 0, 0, 1/2]
r = 2           # number of jumps
max_z = 20      # approximate Gaussian integral from -20σ to 20σ
num_pts = 500   # approximate Gaussian integral with 500 points
>>> parisi_minimize(C_p_squared, r, max_z, num_pts) # outputs Parisi value
```

This code optimizes over piecewise constant functions with $r$ jumps.
Be careful! The runtime is $ \Omega( num_pts^{r+1} ) $.
This example above may take 15 minutes on a laptop.
