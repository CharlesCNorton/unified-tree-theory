import sys
from collections import defaultdict

sys.setrecursionlimit(10**7)

def generate_full_binary_trees(n, memo):
    """
    Return a list of all full binary tree shapes that have n internal nodes.
    Uses memoization to avoid regenerating shapes multiple times.
    """
    if n in memo:
        return memo[n]

    if n == 0:
        return [None]
    if n == 1:
        return [("node", None, None)]

    result = []
    for i in range(n):
        j = (n - 1) - i
        left_shapes = generate_full_binary_trees(i, memo)
        right_shapes = generate_full_binary_trees(j, memo)
        for L in left_shapes:
            for R in right_shapes:
                result.append(("node", L, R))
    memo[n] = result
    return result

def count_cherries(tree):
    if not tree:
        return 0
    _, L, R = tree
    root_cherry = (1 if (L is None and R is None) else 0)
    return root_cherry + count_cherries(L) + count_cherries(R)

def build_cherry_coeff_table(max_n):
    """DP to compute c[n][k], from the corrected logic."""
    c = [defaultdict(int) for _ in range(max_n+1)]
    c[0][0] = 1
    for n in range(1, max_n+1):
        for i in range(n):
            j = (n-1) - i
            for p, valL in c[i].items():
                for q, valR in c[j].items():
                    delta = 1 if (i==0 and j==0) else 0
                    c[n][p + q + delta] += valL * valR
    return c

def test_larger_n(max_n=10):
    """Test enumeration vs DP up to n=10 by default."""
    ctable = build_cherry_coeff_table(max_n)
    memo = {}
    all_good = True

    for n in range(max_n+1):
        shapes = generate_full_binary_trees(n, memo)
        dist = defaultdict(int)
        for shp in shapes:
            k = count_cherries(shp)
            dist[k] += 1

        sum_enum = sum(dist.values())
        sum_dp   = sum(ctable[n].values())
        if sum_enum != sum_dp:
            print(f"[FAIL] n={n}: enumerated total={sum_enum}, DP total={sum_dp}")
            all_good = False

        all_k = set(dist.keys()) | set(ctable[n].keys())
        for k in sorted(all_k):
            e = dist[k]
            d = ctable[n][k]
            if e != d:
                print(f"[FAIL] n={n}, k={k}: enumerated={e}, DP={d}")
                all_good = False

        if all_good:
            print(f"n={n} => OK. (#trees={sum_enum})")
        else:
            print(f"n={n} => see mismatches above.")

    if all_good:
        print(f"\nAll distributions match up to n={max_n}. SUCCESS!")
    else:
        print(f"\nSome mismatches up to n={max_n}. See logs above.")

if __name__ == "__main__":
    test_larger_n(max_n=10)   # Adjust as desired (8, 9, 10, etc.)
