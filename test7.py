import unittest
import sympy

# We use sympy for symbolic computation.
x, u = sympy.symbols('x u')

def ternary_tree_series_coefficients(N, u_sym):
    """
    Compute the first N coefficients (starting at n=1) in the series expansion of
    T(x,u) = x + T(x,u)^3 + (u-1)*x^3,
    where T(x,u) = sum_{n>=1} a_n(u) x^n and a_n(u) are polynomials in u.
    
    The recurrence follows from equating:
      T(x,u) = x + (T(x,u))^3 + (u-1)*x^3.
    Since the term (T(x,u))^3 contributes, for n>=3, a coefficient
      sum_{i+j+k=n} a_i(u)*a_j(u)*a_k(u),
    we can solve recursively for a_n(u).
    """
    # We'll use a list a[1..N] for coefficients; index 0 is unused.
    a = [None] * (N+1)
    a[1] = sympy.Integer(1)  # coefficient for x^1 comes from the leaf: a_1 = 1.
    # For n=2, there is no contribution in the right‐hand side (the cube starts at x^3).
    if N >= 2:
        a[2] = sympy.Integer(0)
    # For n>=3, the equation becomes:
    #   a_n(u) = [coefficient of x^n in T(x,u)^3]  + (u-1) if n == 3, else just the cube term.
    for n in range(3, N+1):
        # Compute the coefficient from T(x,u)^3: sum_{i+j+k=n} a_i * a_j * a_k,
        # where indices i,j,k run over 1..n-? (only terms with i,j,k >= 1 contribute).
        coeff = sympy.Integer(0)
        for i in range(1, n):
            for j in range(1, n - i + 1):
                k = n - i - j
                if k < 1:
                    continue
                coeff += a[i] * a[j] * a[k]
        # If n == 3, we must add the explicit (u-1) contribution.
        if n == 3:
            a[n] = coeff + (u_sym - 1)
        else:
            a[n] = coeff
    return a

class TestTernaryTreeGF(unittest.TestCase):
    def test_series_coefficients_specialization(self):
        """
        Test that our computed series coefficients for T(x,u) specialize correctly.
        In particular, when u = 1 we should recover the generating function for full ternary trees:
          T(x;1) = x + T(x;1)^3.
        The known series for full ternary trees (starting at x^1) is:
          x^1: 1,
          x^2: 0,
          x^3: 1,
          x^4: 0,
          x^5: 3,
          x^6: 0,
          x^7: 12,
          x^8: 0,
          x^9: 55, etc.
        (These are the Fuss–Catalan numbers for m=3.)
        """
        N = 10
        coeffs = ternary_tree_series_coefficients(N, u)
        # Now specialize to u = 1.
        coeffs_u1 = [sympy.simplify(c.subs(u, 1)) for c in coeffs[1:]]  # ignore index 0
        # Expected coefficients: (we include terms up to x^10)
        # a1 = 1, a2 = 0, a3 = 1, a4 = 0, a5 = 3, a6 = 0, a7 = 12, a8 = 0, a9 = 55, a10 = 0
        expected = [1, 0, 1, 0, 3, 0, 12, 0, 55, 0]
        self.assertEqual(coeffs_u1, expected,
                         f"Specialization to u=1 did not match expected series: {coeffs_u1} != {expected}")

    def test_series_coefficients_unified(self):
        """
        Test that for u ≠ 1 the series is different.
        For example, when u = 2, the coefficient a_3 should become: a_3 = (u-1)+1 = u,
        so when u=2 we expect 2.
        """
        N = 5
        coeffs = ternary_tree_series_coefficients(N, u)
        # a_3 should be (u-1)+1 = u, so when u=2 we expect 2.
        a3_u2 = coeffs[3].subs(u, 2)
        # Use nsimplify to get an exact expression:
        a3_u2_simplified = sympy.nsimplify(a3_u2)
        self.assertEqual(a3_u2_simplified, sympy.Integer(2),
                         f"For u=2, expected a_3 = 2, got a_3 = {a3_u2_simplified}")

if __name__ == '__main__':
    unittest.main()
