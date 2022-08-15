# CSP to Parisi

This Python code calculates the average satisfying fraction of a random sparse CSP, using the Parisi formula.

packages required: SymPy, NumPy, SciPy

```bash
git clone git@github.com:marwahaha/csp-parisi.git
cd csp-parisi
```
```bash
python -i csp_parisi.py
```

## Usage

To get the value of a random CSP with Boolean predicate $f$:

```python
k = 3             # number of variables
inp = '01111111'  # specifies 1=T, 0=F for a predicate
get_truth_table(k, inp)
calculate_csp_value(k, inp)
```
When the output is ${0,1}$, the value is the maximum satisfying fraction of the CSP.

Some methods you may also want to use are defined below:


### Get Fourier weights $[|f_{=0}|^2,...,|f_{=k}|^2]$ of a Boolean predicate $f$


```python
k = 3                       # number of variables
inp = '01111111'            # specifies 1=T, 0=F for a predicate
get_truth_table(k, inp)     # to see the truth table
get_fourier_weights(k, inp) # outputs Fourier weights
```

**Note**:
This assumes the predicate has codomain ${-1, +1}$.
If you need ${0,1}$ output:
* Divide all nonzero Fourier weights by 4.
* The weight at zero is the square of the average value of the function. You can use `average_value_zero_one(inp)` to calculate the average value.


### Get approximate Parisi value, given list of mixture function coefficients $[-, c_1^2, c_2^2, ...]$

```bash
python -i get_parisi_value.py
```

```python
C_p_squared = [0, 0, 0, 1/2]
parisi_two_piece_approx(C_p_squared) # outputs Parisi value
```

The mixture function does not include $c_0$, so the 0th entry of `C_p_squared` is ignored.

This code only optimizes over piecewise constant functions with two pieces.
It is reasonably fast (~15 seconds to run), and is >99% accurate on known results.
You can run `verify_parisi_two_piece_accuracy()` to double-check the accuracy claims.


### Get more accurate Parisi value, given list of mixture function coefficients $[-, c_1^2, c_2^2, ...]$

```bash
python -i get_parisi_value.py
```

```python
C_p_squared = [0, 0, 0, 1/2]
r = 2           # number of pieces
max_z = 20      # approximate Gaussian integral from -20σ to 20σ
num_pts = 200   # approximate Gaussian integral with 500 points
parisi_minimize(C_p_squared, r, max_z, num_pts) # outputs Parisi value
```

The mixture function does not include $c_0$, so the 0th entry of `C_p_squared` is ignored.

This code optimizes over piecewise constant functions with $r$ pieces.
Be careful! The runtime is $ \Omega( num_pts^{r+1} ) $.
This example above may take 15 minutes on a laptop.

### Get values for all Boolean predicates

This enumerates through all truth tables.
Running up through $k=3$ takes 5-10 minutes on a laptop.

```python
for k in range(1, 4):
    print("Predicates of", k, "variables")
    for i in range(2**(2**k)):
        inp = get_0101_from_num(2**k, i)
        print(inp)
        calculate_csp_value(k, inp)
```

The following enumerates through symmetric predicates; i.e. those that depend only on Hamming weight.

```python
for k in range(1, 6):
    print("Predicates of", k, "variables")
    lookup = [sum([int(i) for i in get_0101_from_num(k, num)]) for num in range(2**k)]
    for i in range(2**(k+1)):
        on_off_list = get_0101_from_num(k+1, i)
        print('on off:', on_off_list)
        inp = ''.join([on_off_list[i] for i in lookup])
        print('predicate:', inp)
        calculate_csp_value(k, inp)
        print()
```
