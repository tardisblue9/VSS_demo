import random
import math
# from decimal import Decimal
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
    # print("q = " + str(q))
    # print("\nq is prime\n")

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
        # if len(Z_p_star) > 1:
        #     break

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


def generate_shares(N, T, K, secrets, p, q, r, g, alphas, betas):
    secrets_check = []
    for secret in secrets:
        if secret == 0:
            secrets_check.append(1)  # "0": one exception after quantization
        else:
            secrets_check.append(secret)
    secrets = secrets_check
    for secret in secrets:
        assert secret >= 1 and secret <= q, "secret not in range"

    FIELD_SIZE = q
    noises = [random.randrange(0, FIELD_SIZE) for _ in range(T)]

    shares = []
    for alpha in alphas:
        y = _lagrange_interpolate(alpha, betas, secrets + noises, q)
        shares.append((alpha, y))

    commitments = commitment(secrets + noises, g, p)
    start = time.time()

    verifications = []
    for alpha in alphas:
        # check1 = g ** shares[i-1][1] % p
        # check1 = g ** share_ith(shares, i) % p
        check1 = quick_pow(g, share_ith(shares, alpha), p)
        check2 = verification(commitments, alpha, betas, p, q)
        verifications.append(check2)
        if (check1 % p) == (check2 % p):
            pass
        else:
            print("checking fails with:", check1, check2)
            1/0
        print(alpha, "-th user ============= tag at time ", time.time()-start,"seconds =============")
        start = time.time()
    # commitments, verifications = [0,], [0,]
    return shares, commitments, verifications


def share_ith(shares, i):
    for share in shares:
        if share[0] == i:
            return share[1]
    return None


def quick_pow(a, b, q):  # compute a^b mod q, in a faster way？
    temp = 1
    for i in range(1, b + 1):
        temp = (temp * a) % q
    return temp % q


def commitment(paras, g, p):
    commitments = []
    for para in paras:
        # c = g ** coefficient_value % p
        c = quick_pow(g, para, p)
        commitments.append(c)
    return commitments


def verification(commitments, alpha, betas, p, q):
    # v_pos, v_neg = 1, 1
    v = 1
    for i, c in enumerate(commitments):
        num, den = 1, 1
        for k, _ in enumerate(commitments):
            if k != i:
                num *= alpha - betas[k]
                den *= betas[i] - betas[k]
            else:
                pass
        # if num / den > 0:
        #     v_pos = v_pos * quick_pow(c, int(num / den) % q, p) # c ** int(num / den) % p
        # else:
        #     v_neg = v_neg * quick_pow(c, int(- num / den) % q, p) # c ** int(-num / den) % p

        # v = (v * quick_pow(c, int(num / den) % q, p)) % p
        v = (v * quick_pow(c, _divmod(num, den, q) % q , p)) % p
    # v = _divmod(v_pos, v_neg, p)
    return v


def reconstruct_secret(pool, q, betas, K):
    out = []
    x_s, y_s = [], []
    for share in pool:
        x_s.append(int(share[0]))
        y_s.append(int(share[1]))
    for k in range(K):
        beta = betas[k]
        # out.append(f_rec(beta,pool,q))
        out.append(_lagrange_interpolate(beta, x_s, y_s, q))
    return out


def _lagrange_interpolate(x, x_s, y_s, q):
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
    L = 0
    for i in range(k):
        others = list(x_s)
        cur = others.pop(i)
        nums.append(PI(x - o for o in others))
        dens.append(PI(cur - o for o in others))

        L += _divmod(y_s[i] * nums[i], dens[i], q) 
    # den = PI(dens)
    # num = sum([_divmod(nums[i] * den * y_s[i] % q, dens[i], q) for i in range(k)])

    # L = sum( [_divmod(y_s[i] * nums[i], dens[i], q) for i in range(k)] )

    return L % q
    # return _divmod(num, den, q) % q


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
    print("========Main LCC Starts==========")
    p,q,r,g = initial(2* 10**4)

    # Secret taken from the group Z_q* 
    # T, N, K = 2, 10, 3
    T, N, K = 7, 40, 15
    secrets = [random.randint(2,q-1) for i in range(K)]
    print(f'Original Secret: {secrets}')

    # Phase I: Generation of shares
    alphas = list(range(1, 1+N))
    betas = list(range(1+N, N+K+T+1))
    shares, commitments, verifications= generate_shares(N, T, K, secrets, p, q, r, g, alphas, betas)
    print(f'Shares: {", ".join(str(share) for share in shares)}')
    # print(f'Commitments: {", ".join(str(commitment) for commitment in commitments)}')
    # print(f'verifications: {", ".join(str(verification) for verification in verifications)}')

    # Phase II: Secret Reconstruction
    # Picking t shares randomly for reconstruction
    pool = random.sample(shares, T+K)
    print(f'Combining shares: {", ".join(str(share) for share in pool)}')
    secret_reconstructed = reconstruct_secret(pool, q, betas, K)
    print("reconstruct_secret:",secret_reconstructed)
    print(f'Original Secret: {secrets}')

    time_end = time.time()
    print('time cost in second:', time_end-time_start)











