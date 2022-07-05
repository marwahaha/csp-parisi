import numpy as np
from scipy import stats
from scipy.optimize import minimize
from scipy.integrate import quad


def a(qs, l, xiprime):
    """
    a_l^2 = xi'(q_{l+1}) - xi'(q_l)
    """
    return ( xiprime(qs[l+1]) - xiprime(qs[l]) )**0.5


def psi0(qs, ms, xiprime, r):
    """
    Psi_{r+1}(x) = abs(x).
    Psi_l(x) = Log( E[Exp( m_l Psi_{l+1}(x + a_l z) )] )/m_l over z~N(0,1).
    Psi_0(x) = E[Psi_1(x + a_0 z)] over z~N(0,1).

    This uses PDF_INPS as an approximation to a Gaussian integral.
    """
    # Get a_l values.
    a_s = [a(qs, l, xiprime) for l in range(r+1)]

    # Build up grid of Gaussian z_0 to z_r.
    start = sum([a_s[i]*GRID[i] for i in range(len(GRID))])

    # l = r+1
    start = np.abs(start)

    # l = 1...r
    for i in list(range(1, r+1))[::-1]:
        start = np.log(np.sum(np.exp(ms[i]*start)*PDF_INPS, axis=i))/ms[i]

    # l = 0
    start = np.sum(PDF_INPS*start, axis=0)

    return start


def penalty(qs, ms, t_times_xiprimeprime, r):
    """
    0.5 Integral[ f(t) t xi''(t) dt ] from 0 to 1.
    This is equal to 0.5 Sum_{i=0}^k [ m_i Integral[ t xi''(t) dt ] from q_i to q_{i+1} ].
    """
    out = 0
    for i in range(r+1):
        integral = quad(t_times_xiprimeprime, qs[i], qs[i+1])[0]
        out += ms[i] * integral
    return 0.5 * out


def parisi_value(qs, ms, xiprime, t_times_xiprimeprime, r):
    """
    Value of Parisi functional given qs and ms.
    """
    return psi0(qs, ms, xiprime, r) - penalty(qs, ms, t_times_xiprimeprime, r)


def make_parisi_calculator(xiprime, t_times_xiprimeprime, r):
    """
    Creates a function that takes a list of adjustments:
    (m_1, m_2-m_1, ...,m_k-m_{r-1}, q_1, q_2-q_1,...,q_k-q_{r-1})
    and outputs the value of the Parisi functional.

    The function is described as
    0  : [ 0 , q_1)
    m_1: [q_1, q_2)
    m_2: [q_2, q_3)
    ...
    m_r: [q_r, 1)
    """
    def parisi_calculator(inp):
        # assert len(inp) == 2*r
        inp_qs, inp_ms= inp[:r], inp[r:]

        # if bad input, return a large number
        if np.any(np.array(inp) < 0) or sum(inp_qs) > 1:
            return 10000

        qs = np.array([0,*[sum(inp_qs[:i+1]) for i in range(r)],1])
        ms = np.array([0,*[sum(inp_ms[:i+1]) for i in range(r)]])
        output = parisi_value(qs, ms, xiprime, t_times_xiprimeprime, r)
        return output

    return parisi_calculator


def generate_xi_functions(C_p_squared):
    """
    Returns xi, xi', xi'', and lambda t: t*xi''(t).
    """
    xi = lambda x: sum([x**i * C_p_squared[i] for i in range(len(C_p_squared))])
    xiprime = lambda x: sum([i * x**(i-1) * C_p_squared[i] for i in range(1, len(C_p_squared))])
    xiprimeprime = lambda x: sum([i * (i-1) * x**(i-2) * C_p_squared[i] for i in range(2, len(C_p_squared))])
    t_times_xiprimeprime = lambda x: sum([x * i * (i-1) * x**(i-2) * C_p_squared[i] for i in range(2, len(C_p_squared))])
    return xi, xiprime, xiprimeprime, t_times_xiprimeprime


def approximate_gaussian_integral(max_z, num_pts, r):
    """
    Sets precision of Gaussian integral and
    builds globals INPS, PDF_INPS, and GRID.

    If max_z is too small, the convergence is bad,
    especially for higher coefficients.

    Calculating the Parisi value has runtime
    Omega( num_pts^(r+1) ), where r is the number of jumps.
    """

    global INPS, PDF_INPS, GRID

    INPS = np.linspace(-max_z,max_z,num_pts)

    PDF_INPS = stats.norm.pdf(INPS)
    PDF_INPS = PDF_INPS/np.sum(PDF_INPS)
    assert np.allclose(sum(PDF_INPS),1), sum(PDF_INPS)

    GRID = np.meshgrid(*[INPS]*(r+1), indexing='ij')


def parisi_minimize(C_p_squared, r, max_z, num_pts):
    """
    Takes in:
        * coefficients of mixture function, specified as [c_0^2, c_1^2, c_2^2, ...],
        * r (number of jumps)
        * Gaussian integral precision (max_z score and number of points).

    Uses SciPy to minimize the Parisi functional.

    The mixture function should not have a c_0 term,
    so this is automatically zero'd out.

    If the coefficients are too small,
    the convergence is not very good.
    """

    # set Gaussian integral precision
    approximate_gaussian_integral(max_z, num_pts, r)

    # zero out c_0^2
    C_p_squared = C_p_squared.copy()
    C_p_squared[0] = 0

    # make sure input is valid
    for val in C_p_squared:
        assert val >= 0, "all entries of C_p_squared must be non-negative"

    _, xiprime, _, t_times_xiprimeprime = generate_xi_functions(C_p_squared)
    calculator = make_parisi_calculator(xiprime, t_times_xiprimeprime, r)

    # if calculator is zero, then use zero function
    if sum(C_p_squared[:2]) == sum(C_p_squared):
        calculator = lambda *args: 0

    return minimize(calculator,
                    [np.random.random()/r for _ in range(2*r)],
                    method='Powell',
                    options={"xtol": 1e-10, "ftol":1e-14}
                   )


def parisi_one_jump_approx(C_p_squared, print_output = False):
    """
    Finds the minimum Parisi value
    among piecewise constant functions
    with r=1 jump.

    This is a reasonable approximation to known values:
        >99.7% close for Max 2XOR/MaxCut
        >99.95% close for Max 3XOR
    """

    r = 1
    opt = parisi_minimize(C_p_squared, r, 20, 500)

    if print_output:
        print("coeffs:", C_p_squared)
        print("best adjustments:", opt.x)
        print("Parisi value:", opt.fun)

    return float(opt.fun)


def verify_parisi_one_jump_accuracy():
    """
    Code that asserts the accuracy claims made above.

    reference for known values:
    https://arxiv.org/pdf/2009.11481.pdf
    """
    print("Claim: Max 2XOR is 99.7% accurate.")
    assert np.isclose(parisi_one_jump_approx([0, 0, 1/2]), 0.763168, rtol=3e-3), "Max 2XOR is not 99.7% accurate"
    print("Claim: Max 3XOR is 99.95% accurate.")
    assert np.isclose(parisi_one_jump_approx([0, 0, 0, 1/2]), 0.8132, rtol=5e-4), "Max 3XOR is not 99.95% accurate"
    print("Claims are correct.")