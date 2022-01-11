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


def generate_shares(n, t, secret, p, q, r, g, idxs, idx):
    if secret==0: secret=1 # one exception after quantization
    assert secret >= 1 and secret <= q, "secret not in range"

    FIELD_SIZE = q
    coefficients = coeff(t, secret, FIELD_SIZE)
    # users = list(range(1, n+1)) # can be modified
    users = list(idxs)
    users.pop(idx)  # users are to recieve the shares

    shares = []
    for i in users:
        f_i = f(i, coefficients, q)
        shares.append((i, f_i))

    commitments = commitment(coefficients, g, p)
    verifications = []
    for i in users:
        # check1 = g ** shares[i-1][1] % p
        check1 = g ** share_ith(shares, i) % p
        check2 = verification(g, commitments, i, p)
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
    y = 0
    for coefficient_index, coefficient_value in enumerate(coefficients[::-1]):
        y += (x ** coefficient_index * coefficient_value)
    return y


def commitment(coefficients, g, p):
    commitments = []
    for coefficient_index, coefficient_value in enumerate(coefficients[::-1]):
        c = g ** coefficient_value % p
        # c = quick_pow(g,coefficient_value,p)
        commitments.append(c)
    return commitments


def verification(g, commitments, i, p):
    v = 1
    for k, c in enumerate(commitments):
        # v = v * (c) ** (i ** k) % p
        v = v * quick_pow(c,i ** k,p) % p
    return v


def quick_pow(a, b, q):  # compute a^b mod q, in a faster way
    temp = 1
    for i in range(1, b + 1):
        temp = temp * a % q
    return temp % q


def reconstruct_secret(pool, q):
    sums = 0
    prod_arr = []

    for j, share_j in enumerate(pool):
        xj, yj = share_j
        prod = Decimal(1)

        for i, share_i in enumerate(pool):
            xi, _ = share_i
            if i != j:
                prod *= Decimal(xi / (xi - xj))

        prod *= yj
        sums += Decimal(prod)

    out = int(Decimal(sums))
    if out == 1: return 0 # back to the exception mentioned above
    return out
    # return int(round(Decimal(sums)))


# Driver code
if __name__ == '__main__':
    # initialization
    time_start = time.time()
    print("========Main VSS Starts==========")
    p,q,r,g = initial(10**3)

    # Secret taken from the group Z_q* 
    t, n = 2, 9
    secrets = [47,83,77,89,99,78,89,90,13]
    for secret in secrets:
        print(f'Original Secret: {secret}')
    print('directly added sum:',sum(secrets))

    # Phase I: Generation of shares
    users = [9, 3, 7, 8, 4, 6, 5, 2, 1, 10]
    user = users.index(5)

    shares_list = []
    for secret in secrets:
        shares, commitments, verifications= generate_shares(n, t, secret, p, q, r, g, users, user)
        shares_list.append(shares)

    # aggregation
    users.pop(user)
    idx_pool = random.sample(users, t)
    pool =[]
    for idx in idx_pool:
        add = 0
        for shares in shares_list:
            for share in shares:
                if share[0] == idx:
                    add += share[1]
                else: pass
        pool.append((idx, add))


    # Phase II: Secret Reconstruction
    # Picking t shares randomly for reconstruction
    print(f'Combining shares: {", ".join(str(share) for share in pool)}')
    secret_reconstructed = reconstruct_secret(pool, q)
    print("reconstruct_secret:",secret_reconstructed)

    time_end = time.time()
    print('time cost in second:', time_end-time_start)












