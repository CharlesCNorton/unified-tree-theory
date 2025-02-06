"""
Comprehensive Test Suite for the Unified Generating Function Framework

This script verifies the following for full binary trees with n leaves:
  - T(n): number of full binary trees (Catalan numbers)
  - S(n): total Sackin index (sum of leaf depths)
  - C(n): total Colless index
  - Phi(n): total cophenetic index
  - X(n): total cherry count
  - S2(n): total Sackin2 index (sum of squares of leaf depths)

It also tests:
  - That the DP recurrences agree with known closed-form formulas (or known small-n values)
  - That the DP values for S(n) match those computed by brute-force enumeration (for small n)
  - That the generating function for the Catalan numbers, when expanded symbolically, yields the same coefficients
  - That the asymptotic behavior of T(n) agrees with the classical estimate
  - That the runtime scales quadratically (via log–log regression)

"""

import math
import time
import unittest
import statistics
import numpy as np
from sympy import symbols, series, sqrt

# ------------------------------
# 1. Dynamic Programming Routines for T, S, C, Phi, X
# ------------------------------

def compute_invariants(N):
    """
    Compute the following for full binary trees with n leaves, for 1 <= n <= N:
      T[n]: number of full binary trees (Catalan numbers)
      S[n]: total Sackin index (sum of leaf depths)
      C[n]: total Colless index
      Phi[n]: total cophenetic index
      X[n]: total cherry count
    Returns five lists T, S, C, Phi, X indexed from 1 to N.
    """
    T = [0] * (N + 1)      # Catalan numbers
    S = [0] * (N + 1)      # Sackin index
    C = [0] * (N + 1)      # Colless index
    Phi = [0] * (N + 1)    # Total cophenetic index
    X = [0] * (N + 1)      # Cherry count

    # Base cases
    T[1] = 1
    S[1] = 0
    C[1] = 0
    Phi[1] = 0
    X[1] = 0

    if N >= 2:
        T[2] = 1       # one tree with 2 leaves
        S[2] = 2       # both leaves at depth 1: 1+1=2
        C[2] = 0       # imbalance: |1-1|=0
        Phi[2] = 0     # only pair: LCA at depth 0 -> 0
        X[2] = 1       # the unique tree is a cherry

    # Helper function: binom(n, 2)
    def binom2(n):
        return n * (n - 1) // 2

    for n in range(3, N + 1):
        T_n = 0
        S_n = 0
        C_n = 0
        Phi_n = 0
        X_n = 0
        for i in range(1, n):
            j = n - i
            prod = T[i] * T[j]
            T_n += prod
            S_n += S[i] * T[j] + S[j] * T[i] + n * prod
            C_n += C[i] * T[j] + C[j] * T[i] + abs(2 * i - n) * prod
            Phi_n += Phi[i] * T[j] + Phi[j] * T[i] + (binom2(i) + binom2(j)) * prod
            X_n += X[i] * T[j] + X[j] * T[i]
        T[n] = T_n
        S[n] = S_n
        C[n] = C_n
        Phi[n] = Phi_n
        X[n] = X_n

    return T, S, C, Phi, X

def compute_T(n): return compute_invariants(n)[0][n]
def compute_S(n): return compute_invariants(n)[1][n]
def compute_C(n): return compute_invariants(n)[2][n]
def compute_Phi(n): return compute_invariants(n)[3][n]
def compute_X(n): return compute_invariants(n)[4][n]

# ------------------------------
# 2. Dynamic Programming Routine for S2 (Sackin2 Index)
# ------------------------------
def compute_invariants_S2(N):
    """
    Compute T[n] and S2[n] for full binary trees with n leaves.
    S2(n) = sum_{T with n leaves} (sum of squares of leaf depths).
    Recurrence:
      S2(1)=0.
      For n>=2:
        S2(n)= sum_{i=1}^{n-1} [ S2(i)*T(n-i) + S2(n-i)*T(i)
                  + 2*(S(i)*T(n-i) + S(n-i)*T(i)) + n*T(i)*T(n-i) ]
    """
    T = [0] * (N + 1)
    S = [0] * (N + 1)
    S2 = [0] * (N + 1)

    # Base cases
    T[1] = 1
    S[1] = 0
    S2[1] = 0

    if N >= 2:
        T[2] = 1
        S[2] = 2   # as before
        S2[2] = 2  # unique tree with two leaves: 1^2+1^2=2

    for n in range(3, N + 1):
        T_n = 0
        S_n = 0
        S2_n = 0
        for i in range(1, n):
            j = n - i
            prod = T[i] * T[j]
            T_n += prod
            S_n += S[i] * T[j] + S[j] * T[i] + n * prod
            S2_n += (S2[i] * T[j] + S2[j] * T[i] +
                     2 * (S[i] * T[j] + S[j] * T[i]) +
                     n * prod)
        T[n] = T_n
        S[n] = S_n
        S2[n] = S2_n

    return T, S, S2

def compute_S2(n): return compute_invariants_S2(n)[2][n]

# ------------------------------
# 3. Brute-Force Generation and Sackin Index Computation (for small n)
# ------------------------------
def generate_trees(n):
    """
    Generate all full binary trees with n leaves.
    Represent a leaf as 'L' and an internal node as a tuple (left, right).
    """
    if n == 1:
        return ['L']
    trees = []
    for i in range(1, n):
        j = n - i
        left_trees = generate_trees(i)
        right_trees = generate_trees(j)
        for left in left_trees:
            for right in right_trees:
                trees.append((left, right))
    return trees

def compute_tree_S(tree, depth=0):
    """
    Compute the Sackin index of a tree (sum of leaf depths).
    """
    if tree == 'L':
        return depth
    left, right = tree
    return compute_tree_S(left, depth + 1) + compute_tree_S(right, depth + 1)

def brute_force_compute_S(n):
    trees = generate_trees(n)
    total = 0
    for t in trees:
        total += compute_tree_S(t)
    return total

# ------------------------------
# 4. Unit Testing with unittest
# ------------------------------
class TestAllInvariants(unittest.TestCase):

    def test_catalan_numbers(self):
        # Test T(n) against the closed-form Catalan numbers: T(n) = C_{n-1} = (1/n)*comb(2n-2, n-1)
        for n in range(1, 21):
            T_computed = compute_T(n)
            T_closed = math.comb(2*(n-1), n-1) // n
            self.assertEqual(T_computed, T_closed)

    def test_sackin_index(self):
        # Known small values: S(1)=0, S(2)=2, S(3)=10, S(4)=44, S(5)=186
        expected = {1: 0, 2: 2, 3: 10, 4: 44, 5: 186}
        for n, exp_val in expected.items():
            self.assertEqual(compute_S(n), exp_val)

    def test_colless_index(self):
        # Known small values: C(1)=0, C(2)=0, C(3)=2, C(4)=12, C(5)=62
        expected = {1: 0, 2: 0, 3: 2, 4: 12, 5: 62}
        for n, exp_val in expected.items():
            self.assertEqual(compute_C(n), exp_val)

    def test_cophenetic_index(self):
        # Known small values: Phi(1)=0, Phi(2)=0, Phi(3)=2, Phi(4)=18, Phi(5)=116
        expected = {1: 0, 2: 0, 3: 2, 4: 18, 5: 116}
        for n, exp_val in expected.items():
            self.assertEqual(compute_Phi(n), exp_val)

    def test_cherry_count(self):
        # Known small values: X(1)=0, X(2)=1, X(3)=2, X(4)=6, X(5)=20
        expected = {1: 0, 2: 1, 3: 2, 4: 6, 5: 20}
        for n, exp_val in expected.items():
            self.assertEqual(compute_X(n), exp_val)

    def test_sackin2_index(self):
        # Expected values computed manually:
        # n=1: 0; n=2: 2; n=3: 18; n=4: 108; n=5: 562
        expected = {1: 0, 2: 2, 3: 18, 4: 108, 5: 562}
        for n, exp_val in expected.items():
            self.assertEqual(compute_S2(n), exp_val)

    def test_brute_force_vs_dp_sackin(self):
        # For small n (n <= 6), compare brute-force Sackin index with DP result.
        for n in range(1, 7):
            dp_val = compute_S(n)
            bf_val = brute_force_compute_S(n)
            self.assertEqual(dp_val, bf_val)
            trees = generate_trees(n)
            self.assertEqual(len(trees), compute_T(n))

    def test_symbolic_series_catalan(self):
        # Expand the generating function G(x) = (1 - sqrt(1-4*x)) / 2 and compare coefficients.
        x = symbols('x')
        G = (1 - sqrt(1 - 4*x)) / 2
        G_series = series(G, x, 0, 12).removeO()
        for n in range(1, 12):
            coeff_sym = G_series.coeff(x, n)
            coeff = int(coeff_sym)
            self.assertEqual(compute_T(n), coeff)

    def test_asymptotic_catalan(self):
        # Check that T(n) ~ 4^n/(4*sqrt(pi)*n^(3/2)) for n=50,100 (within 10% relative error).
        for n in [50, 100]:
            T_n = compute_T(n)
            asymptotic = 4**n / (4 * math.sqrt(math.pi) * n**(1.5))
            relative_error = abs(T_n - asymptotic) / asymptotic
            self.assertLess(relative_error, 0.10)

    def test_runtime_scaling(self):
        # Measure runtime for computing S(n) for various n, average over 5 runs, and perform log-log regression.
        sample_ns = [50, 75, 100, 125, 150]
        avg_times = []
        num_runs = 5
        for n in sample_ns:
            times = []
            for _ in range(num_runs):
                t0 = time.time()
                compute_S(n)
                times.append(time.time() - t0)
            avg_times.append(statistics.mean(times))
        log_ns = np.log(sample_ns)
        log_times = np.log(avg_times)
        slope, intercept = np.polyfit(log_ns, log_times, 1)
        self.assertTrue(1.8 < slope < 2.2, msg=f"Expected quadratic scaling (slope ~2) but got slope = {slope:.2f}")

# ------------------------------
# 5. Main Block to Run Tests
# ------------------------------
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestAllInvariants)
    runner = unittest.TextTestRunner()
    runner.run(suite)
