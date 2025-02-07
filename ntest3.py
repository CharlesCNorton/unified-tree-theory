import math
import sympy as sp
import pandas as pd

# Disable logging for raw output.
import logging
logging.basicConfig(level=logging.CRITICAL)

# ---------------------------
# Recurrence-based computation for invariants.
# ---------------------------
def compute_invariants(n_max):
    T = [0] * (n_max + 1)
    S = [0] * (n_max + 1)
    C = [0] * (n_max + 1)
    Phi = [0] * (n_max + 1)
    X = [0] * (n_max + 1)
    S2 = [0] * (n_max + 1)

    T[1] = 1
    S[1] = 0
    C[1] = 0
    Phi[1] = 0
    X[1] = 0
    S2[1] = 0

    if n_max >= 2:
        T[2] = 1
        S[2] = 2
        C[2] = 0
        Phi[2] = 0
        X[2] = 1
        S2[2] = 2

    for n in range(3, n_max + 1):
        T_n = 0
        S_n = 0
        C_n = 0
        Phi_n = 0
        X_n = 0
        S2_n = 0
        for i in range(1, n):
            j = n - i
            prod = T[i] * T[j]
            T_n += prod
            S_n += S[i] * T[j] + S[j] * T[i] + n * prod
            C_n += C[i] * T[j] + C[j] * T[i] + abs(2 * i - n) * prod
            binom_i = i * (i - 1) // 2
            binom_j = j * (j - 1) // 2
            Phi_n += Phi[i] * T[j] + Phi[j] * T[i] + (binom_i + binom_j) * prod
            X_n += X[i] * T[j] + X[j] * T[i]
            S2_n += S2[i] * T[j] + S2[j] * T[i] + 2 * (S[i] * T[j] + S[j] * T[i]) + n * prod
        T[n] = T_n
        S[n] = S_n
        C[n] = C_n
        Phi[n] = Phi_n
        X[n] = X_n
        S2[n] = S2_n
    return T, S, C, Phi, X, S2

# ---------------------------
# Closed-form generating function series expansion.
# ---------------------------
def series_coefficients(expr, n_max):
    x = sp.symbols('x')
    poly = sp.series(expr, x, 0, n_max + 1).removeO()
    coeffs = [sp.expand(poly).coeff(x, n) for n in range(n_max + 1)]
    return coeffs

# ---------------------------
# Catalan number via formula.
# ---------------------------
def catalan(n):
    return math.comb(2*n, n) // (n + 1)

# ---------------------------
# Main testing function.
# ---------------------------
def test_all(n_max=20):
    x = sp.symbols('x')

    # Closed-form generating functions:
    T_closed_expr = (1 - sp.sqrt(1 - 4*x)) / 2
    Q_closed_expr = x * (1 - sp.sqrt(1 - 4*x)) / (1 - 4*x)  # Sackin index generating function.
    P_closed_expr = x * ((1 - 4*x)**(sp.Rational(3, 2)) - 1 + 6*x - 4*x**2) / (2 * (1 - 4*x)**(sp.Rational(3, 2)))  # Corrected Colless.
    U_closed_expr = (4*x*(1 - sp.sqrt(1 - 4*x) - 2*x)) / ((1 - 4*x)**(sp.Rational(3, 2))) + (x*(1 - sp.sqrt(1 - 4*x))) / (1 - 4*x)  # Sackin2 index.

    # Series expansion for closed-form:
    T_series = series_coefficients(T_closed_expr, n_max)
    Q_series = series_coefficients(Q_closed_expr, n_max)
    P_series = series_coefficients(P_closed_expr, n_max)
    U_series = series_coefficients(U_closed_expr, n_max)

    # Recurrence-based values.
    T_rec, S_rec, C_rec, Phi_rec, X_rec, S2_rec = compute_invariants(n_max)

    # Table 1: Recurrence vs Closed Form for T, S, C, S2.
    data1 = []
    for n in range(1, n_max + 1):
        row = {
            "n": n,
            "T(n) rec": T_rec[n],
            "T(n) closed": int(T_series[n]),
            "S(n) rec": S_rec[n],
            "S(n) closed": int(Q_series[n]),
            "C(n) rec": C_rec[n],
            "C(n) closed": int(P_series[n]),
            "S2(n) rec": S2_rec[n],
            "S2(n) closed": int(U_series[n])
        }
        data1.append(row)
    df1 = pd.DataFrame(data1)
    
    # Table 2: Catalan check: Compare T(n) rec with catalan(n-1).
    data2 = []
    for n in range(1, n_max + 1):
        if n == 1:
            cat = 1
        else:
            cat = catalan(n - 1)
        data2.append({"n": n, "T(n) rec": T_rec[n], "Catalan(n-1)": cat, "Difference": T_rec[n] - cat})
    df2 = pd.DataFrame(data2)
    
    # Table 3: Averages: E[S]=S(n)/T(n), E[C]=C(n)/T(n), E[S2]=S2(n)/T(n), E[X]=X(n)/T(n)
    data3 = []
    for n in range(1, n_max + 1):
        avg_S = S_rec[n] / T_rec[n] if T_rec[n] != 0 else 0
        avg_C = C_rec[n] / T_rec[n] if T_rec[n] != 0 else 0
        avg_S2 = S2_rec[n] / T_rec[n] if T_rec[n] != 0 else 0
        avg_X = X_rec[n] / T_rec[n] if T_rec[n] != 0 else 0
        data3.append({"n": n, "E[S]": avg_S, "E[C]": avg_C, "E[S2]": avg_S2, "E[X]": avg_X})
    df3 = pd.DataFrame(data3)
    
    # Table 4: Asymptotic test for T(n): T(n) ~ 4^n/(4√(π)n^(3/2))
    data4 = []
    for n in range(10, n_max + 1, 5):
        T_val = T_rec[n]
        asym_T = 4**n / (4 * math.sqrt(math.pi) * n**1.5)
        ratio = T_val / asym_T
        data4.append({"n": n, "T(n) rec": T_val, "Asymptotic T(n)": asym_T, "Ratio": ratio})
    df4 = pd.DataFrame(data4)

    # Table 5: Phi(n) (only recurrence values)
    data5 = []
    for n in range(1, n_max + 1):
        data5.append({"n": n, "Phi(n) rec": Phi_rec[n]})
    df5 = pd.DataFrame(data5)

    print("Table 1: Recurrence vs Closed Form (T, S, C, S2)")
    print(df1.to_string(index=False))
    print("\nTable 2: Catalan Check (T(n) rec vs Catalan(n-1))")
    print(df2.to_string(index=False))
    print("\nTable 3: Average Invariants (E[S], E[C], E[S2], E[X])")
    print(df3.to_string(index=False))
    print("\nTable 4: Asymptotic Test for T(n)")
    print(df4.to_string(index=False))
    print("\nTable 5: Total Cophenetic Index Phi(n) (Recurrence-based)")
    print(df5.to_string(index=False))

if __name__ == "__main__":
    test_all(20)
