import math

# Set maximum number of leaves
N = 100

# Initialize arrays (1-indexed: index n corresponds to trees with n leaves)
T = [0] * (N + 1)      # T[n]: number of full binary trees with n leaves
S = [0] * (N + 1)      # S[n]: total Sackin index
C = [0] * (N + 1)      # C[n]: total Colless index
Phi = [0] * (N + 1)    # Phi[n]: total cophenetic index
X = [0] * (N + 1)      # X[n]: total number of cherries

# Base cases
T[1] = 1
S[1] = 0
C[1] = 0
Phi[1] = 0
X[1] = 0

# For n = 2 (the unique tree with two leaves)
if N >= 2:
    T[2] = 1      # Catalan(1) = 1
    S[2] = 2      # Both leaves at depth 1: 1 + 1 = 2
    C[2] = 0      # Imbalance: |1 - 1| = 0
    Phi[2] = 0    # Only one pair of leaves, LCA at depth 0: 0
    X[2] = 1      # The root is a cherry

# Define a helper function for binomial(n,2)
def binom2(n):
    return n * (n - 1) // 2

# Compute recurrences for n from 3 to N
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
        
        # Sackin index: add contributions from left and right, plus extra depth n for each leaf
        S_n += S[i] * T[j] + S[j] * T[i] + n * prod
        
        # Colless index: extra imbalance at root is |2*i - n|
        C_n += C[i] * T[j] + C[j] * T[i] + abs(2 * i - n) * prod
        
        # Cophenetic index: extra contribution at root is binom2(i) + binom2(j)
        Phi_n += Phi[i] * T[j] + Phi[j] * T[i] + (binom2(i) + binom2(j)) * prod
        
        # Cherry count: for n>=3, the extra cherry occurs only when i == j == 1 (which would imply n=2)
        X_n += X[i] * T[j] + X[j] * T[i]
    T[n] = T_n
    S[n] = S_n
    C[n] = C_n
    Phi[n] = Phi_n
    X[n] = X_n

# Print a table of results for n from 1 to N
print("n\tT(n)\t\tS(n)\t\tC(n)\t\tPhi(n)\t\tX(n)")
for n in range(1, N + 1):
    # Use scientific notation for very long integers if needed
    T_str = str(T[n]) if len(str(T[n])) < 35 else f"{T[n]:.5e}"
    S_str = str(S[n]) if len(str(S[n])) < 35 else f"{S[n]:.5e}"
    C_str = str(C[n]) if len(str(C[n])) < 35 else f"{C[n]:.5e}"
    Phi_str = str(Phi[n]) if len(str(Phi[n])) < 35 else f"{Phi[n]:.5e}"
    X_str = str(X[n])
    print(f"{n}\t{T_str}\t{S_str}\t{C_str}\t{Phi_str}\t{X_str}")
