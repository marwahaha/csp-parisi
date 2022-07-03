import sympy

def get_0101_from_num(k, num):
    """
    Convert from a number to a
    binary string of 1=T, 0=F.
    """
    assert 0 <= num, "num must be at least zero"
    assert num < 2**k, "num must be at most 2^k"
    return ('0'*k + bin(num)[2:])[-k:]

assert get_0101_from_num(4, 10) == '1010'


def get_truth_table(k, inp):
    """
    Get truth table from
    the input value.
    """
    assert len(inp) == 2**k, "inp must be exactly 2^k length"
    return {get_0101_from_num(k, i): int(inp[i]) for i in range(len(inp))}

assert get_truth_table(4, '0'*15 + '1') == { '0000': 0,
                                    '0001': 0,
                                    '0010': 0,
                                    '0011': 0,
                                    '0100': 0,
                                    '0101': 0,
                                    '0110': 0,
                                    '0111': 0,
                                    '1000': 0,
                                    '1001': 0,
                                    '1010': 0,
                                    '1011': 0,
                                    '1100': 0,
                                    '1101': 0,
                                    '1110': 0,
                                    '1111': 1 }


def generate_variables(k):
    """
    Generate Sympy variables, from x1 to xk.
    """
    return [sympy.Symbol('x' + str(i+1)) for i in range(k)]


def get_indicator_poly(k, vals, xs):
    """
    Gets indicator polynomial given
    vals (specified as binary string).
    """
    assert len(vals) == k, "vals must be exactly k length"
    assert len(xs) == k, "xs must have exactly k variables"
    out = 1
    for i in range(len(vals)):
        out *= (1 + (2*int(vals[i]) - 1) * xs[i])
    return out/2**k

assert str(get_indicator_poly(4, '1001', generate_variables(4))) == '(1 - x2)*(1 - x3)*(x1 + 1)*(x4 + 1)/16'


def get_fourier_poly(k, inp):
    """
    Takes inp and gets the polynomial
    representing the truth table of inp,
    combining the indicator polynomials.
    Each coefficient is hat{f} on that subset.
    """
    assert len(inp) == 2**k, "inp must be exactly 2^k length"
    xs = generate_variables(k)
    out = 0
    for i in range(len(inp)):
        out += (2*int(inp[i]) - 1) * get_indicator_poly(k, get_0101_from_num(k, i), xs)
    return sympy.expand(sympy.simplify(sympy.expand(out)))

assert str(get_fourier_poly(3, '00100000')) == 'x1*x2*x3/4 - x1*x2/4 + x1*x3/4 - x1/4 - x2*x3/4 + x2/4 - x3/4 - 3/4'


def get_fourier_levels_from_poly(k, fourier_poly):
    """
    Returns |f_{=j}|^2 for each j,
    for 0 <= j <= k,
    given the Fourier polynomial of f.
    """
    if fourier_poly.is_number:
        return [1] + [0,]*k
    poly = sympy.Poly(fourier_poly)
    out = [0]*(k+1)
    for (m, c) in zip(poly.monoms(), poly.coeffs()):
        out[sum(m)] += c**2
    assert sum(out) == 1, "fourier levels must sum to 1"
    return out

assert get_fourier_levels_from_poly(4, get_fourier_poly(4, '1111111111111110')) == [49/64, 1/16, 3/32, 1/16, 1/64]


def get_fourier_levels(k, inp):
    return get_fourier_levels_from_poly(k, get_fourier_poly(k, inp))
