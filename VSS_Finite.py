import random
import math
from decimal import Decimal
import time

def isprime(n):
    if n == 2:
        return True
    if n == 1 or n % 2 == 0:
        return False
    i = 3
    while i <= math.sqrt(n):
        if n % i == 0:
            return False
        i = i + 2
    return True


def initial(Z_lower=100):
    # generate q bigger than z_lower
    q = Z_lower
    while True:
        if isprime(q):
            break
        else:
            q = q + 1
    print("q = " + str(q))
    print("\nq is prime\n")

    # Find p and r
    r = 1
    while True:
        p = r * q + 1
        if isprime(p):
            print("r = " + str(r))
            print("p = " + str(p))
            print("\np is prime\n")
            break
        r = r + 1

    # Compute elements of Z_p*
    Z_p_star = []
    for i in range(0, p):
        if (math.gcd(i, p) == 1):
            Z_p_star.append(i)

    # print("Z_p* = ")
    # print(Z_p_star) # , len(Z_p_star) same length, i.e. range(p)

    # Compute elements of G = {h^r mod p | h in Z_p*}
    G = []
    for i in Z_p_star:
        G.append(i ** r % p)

    G = list(set(G))
    G.sort()
    # print("\nG = ")
    # print(G)
    # print("Order of G is " + str(len(G)) + ". This must be equal to q.")

    # Since the order of G is prime, any element of G except 1 is a generator
    g = random.choice(list(filter(lambda g: g != 1, G)))
    print("\ng = " + str(g) + "\n")

    return p, q, r, g

def generate_shares(n, t, secret, p, q, r, g, idxs):
    if secret==0: secret=1 # one exception after quantization
    assert secret >= 1 and secret <= q, "secret not in range"

    FIELD_SIZE = q
    coefficients = coeff(t, secret, FIELD_SIZE)
    users = list(idxs) # users are to recieve the shares
    assert n==len(users), "these two number should be identical"

    shares = []
    for i in users:
        f_i = f(i, coefficients, q)
        shares.append((i, f_i))

    commitments = commitment(coefficients, g, p)
    verifications = []
    for i in users:
        # check1 = g ** shares[i-1][1] % p
        # check1 = g ** share_ith(shares, i) % p
        check1 = quick_pow(g, share_ith(shares, i), p)
        check2 = verification(g, commitments, int(i), p)
        verifications.append(check2)
        if check1 == check2:
            pass
        else:
            print("checking fails with:", check1, check2)

    return shares, commitments, verifications

def share_ith(shares, i):
    for share in shares:
        if share[0] == i:
            return share[1]
    return None

def coeff(t, secret, FIELD_SIZE):
    coeff = [random.randrange(0, FIELD_SIZE) for _ in range(t - 1)]
    coeff.append(secret)  # a0 is secret
    return coeff

def f(x, coefficients, q):
    y = Decimal('0')
    for coefficient_index, coefficient_value in enumerate(coefficients[::-1]):
        y += (Decimal(str(x)) ** Decimal(str(coefficient_index)) * Decimal(str(coefficient_value)))
        y = Decimal(int(y)%q)
    return int(y)


def commitment(coefficients, g, p):
    commitments = []
    for coefficient_index, coefficient_value in enumerate(coefficients[::-1]):
        # c = g ** coefficient_value % p
        c = quick_pow(g,coefficient_value,p)
        commitments.append(c)
    return commitments


def verification(commitments, i, p):
    v = 1
    for k, c in enumerate(commitments):
        v = v * (c) ** (i ** k) % p
        # v = v * quick_pow(c,i ** k,p) % p
    return v


def quick_pow(a, b, q):  # compute a^b mod q, in a faster way
    temp = 1
    for i in range(1, b + 1):
        temp = temp * a % q
    return temp % q


def reconstruct_secret(pool, q):
    x_s,y_s = [],[]
    for share in pool:
        x_s.append(share[0])
        y_s.append(share[1])
    out = _lagrange_interpolate(0, x_s, y_s, q)
    return out

def _lagrange_interpolate(x, x_s, y_s, p):
    """
    Find the y-value for the given x, given n (x, y) points;
    k points will define a polynomial of up to kth order.
    """
    k = len(x_s)
    assert k == len(set(x_s)), "points must be distinct"
    def PI(vals):  # upper-case PI -- product of inputs
        accum = 1
        for v in vals:
            accum *= v
        return accum
    nums = []  # avoid inexact division
    dens = []
    for i in range(k):
        others = list(x_s)
        cur = others.pop(i)
        nums.append(PI(x - o for o in others))
        dens.append(PI(cur - o for o in others))
    den = PI(dens)
    num = sum([_divmod(nums[i] * den * y_s[i] % p, dens[i], p)
               for i in range(k)])
    return (_divmod(num, den, p) + p) % p

def _extended_gcd(a, b):
    """
    Division in integers modulus p means finding the inverse of the
    denominator modulo p and then multiplying the numerator by this
    inverse (Note: inverse of A is B such that A*B % p == 1) this can
    be computed via extended Euclidean algorithm
    http://en.wikipedia.org/wiki/Modular_multiplicative_inverse#Computation
    """
    x = 0
    last_x = 1
    y = 1
    last_y = 0
    while b != 0:
        quot = a // b
        a, b = b, a % b
        x, last_x = last_x - quot * x, x
        y, last_y = last_y - quot * y, y
    return last_x, last_y

def _divmod(num, den, p):
    """Compute num / den modulo prime p

    To explain what this means, the return value will be such that
    the following is true: den * _divmod(num, den, p) % p == num
    """
    inv, _ = _extended_gcd(den, p)
    return num * inv

# Driver code
if __name__ == '__main__':
    # initialization
    time_start = time.time()
    print("========Main VSS Starts==========")
    p,q,r,g = initial(10**5)

    # Secret taken from the group Z_q* 
    t, n = 5, 10
    secret = 1243
    print(f'Original Secret: {secret}')

    # Phase I: Generation of shares
    users = [9, 3, 7, 8, 4, 6, 5, 2, 1, 10]
    print(n, t, secret, p, q, r, g, users)
    shares, commitments, verifications= generate_shares(n, t, secret, p, q, r, g, users)
    print(f'Shares: {", ".join(str(share) for share in shares)}')
    print(f'Commitments: {", ".join(str(commitment) for commitment in commitments)}')
    print(f'verifications: {", ".join(str(verification) for verification in verifications)}')

    # Phase II: Secret Reconstruction
    # Picking t shares randomly for reconstruction
    pool = random.sample(shares, t)
    print(f'Combining shares: {", ".join(str(share) for share in pool)}')
    secret_reconstructed = reconstruct_secret(pool, q)
    print("reconstruct_secret:",secret_reconstructed)

    time_end = time.time()
    print('time cost in second:', time_end-time_start)












