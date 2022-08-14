from get_fourier_weights import *
from get_parisi_value import *


def average_value_pm_one(inp):
    """
    Gets average value of Boolean predicate,
    assuming output is -1 or +1.
    """
    return sum([2*int(i) - 1 for i in inp])/len(inp)


def average_value_zero_one(inp):
    """
    Gets average value of Boolean predicate,
    assuming output is 0 or 1.
    """
    return sum([int(i) for i in inp])/len(inp)


def calculate_csp_value(k, inp):
    """
    Gets the maximum score (the "value")
    of a CSP constructed by a Boolean predicate,
    using an approximate Parisi value calculation.
    """

    weights = get_fourier_weights(k, inp)

    # The mixture function of the spin glass model
    # is determined by the Fourier weights of the CSP.
    value = parisi_one_piece_approx(weights)

    avg_pm_one = average_value_pm_one(inp)
    print("For -1/+1 output, with a*n clauses:",
            avg_pm_one, "+", value, "/sqrt(a)")

    avg_zero_one = average_value_zero_one(inp)
    print("For  0/1  output, with a*n clauses:",
            avg_zero_one, "+", value/2, "/sqrt(a)")

    return value
