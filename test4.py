"""
Unified Framework Test Suite for Tree Invariants:
Joint Generating Functions, Multivariate Asymptotics, and Extensions

Authors: Charles C. Norton & OpenAI's o3-mini-high
Date: February 5th, 2025

This test suite verifies our unified framework by:
1. Testing full binary tree generation and invariant computations (cherries, Colless, Sackin).
2. Constructing the joint generating function GF_joint(x, y, z) for full binary trees.
3. Specializing the joint GF (setting y = 1, z = 1) and comparing its series expansion,
   truncated to the same order, with the classical Catalan generating function.
4. Running asymptotic tests for F(x,1) and comparing exact coefficients with asymptotic estimates.
5. Testing an extension to full ternary trees.
6. Printing detailed outputs for independent verification.
"""

import sympy as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =============================================================================
# Part I. Utility Functions: Full Binary Trees & Invariants
# =============================================================================

def generate_full_binary_trees(n):
    """Recursively generate all full binary trees with n leaves.
    A leaf is "L"; an internal node is a tuple (left, right)."""
    if n == 1:
        return ["L"]
    trees = []
    for i in range(1, n):
        left_subtrees = generate_full_binary_trees(i)
        right_subtrees = generate_full_binary_trees(n - i)
        for left in left_subtrees:
            for right in right_subtrees:
                trees.append((left, right))
    return trees

def catalan_number(m):
    """Return the m-th Catalan number."""
    return sp.binomial(2*m, m) // (m+1)

def count_leaves(tree):
    """Return the number of leaves in a tree."""
    if tree == "L":
        return 1
    left, right = tree
    return count_leaves(left) + count_leaves(right)

def count_cherries(tree):
    """Return the number of cherries in a full binary tree."""
    if tree == "L":
        return 0
    left, right = tree
    is_cherry = 1 if (left == "L" and right == "L") else 0
    return is_cherry + count_cherries(left) + count_cherries(right)

def colless_index(tree):
    """Compute the Colless index of a full binary tree."""
    if tree == "L":
        return 0
    left, right = tree
    L = count_leaves(left)
    R = count_leaves(right)
    return abs(L - R) + colless_index(left) + colless_index(right)

def sackin_index(tree, depth=0):
    """Compute the Sackin index of a full binary tree (sum of leaf depths)."""
    if tree == "L":
        return depth
    left, right = tree
    return sackin_index(left, depth+1) + sackin_index(right, depth+1)

# =============================================================================
# Part II. Joint Generating Function for Full Binary Trees
# =============================================================================

def joint_generating_function(n, tree_generator):
    """
    For a given n, compute the joint generating function contribution:
      GF = sum_{T in Trees(n)} x^(# leaves) * y^(# cherries) * z^(Colless index).
    """
    x, y, z = sp.symbols('x y z')
    GF = 0
    trees = tree_generator(n)
    for tree in trees:
        GF += x**(count_leaves(tree)) * y**(count_cherries(tree)) * z**(colless_index(tree))
    return sp.simplify(GF)

def aggregate_joint_GF(max_n, tree_generator):
    """
    Aggregate the joint generating function for n = 1 to max_n.
    Returns a sympy expression in variables x, y, z.
    """
    x, y, z = sp.symbols('x y z')
    total_GF = 0
    for n in range(1, max_n+1):
        total_GF += joint_generating_function(n, tree_generator)
    return sp.simplify(total_GF)

# =============================================================================
# Part III. Asymptotic Analysis for Univariate Generating Function
# =============================================================================

def test_asymptotic_estimates():
    """Test that F(x,1) = (1-sqrt(1-4*x))/2 has coefficients matching known asymptotics."""
    x = sp.symbols('x')
    F_univ = (1 - sp.sqrt(1-4*x))/2
    print("\nTesting asymptotic estimates for F(x,1):")
    for n in range(4, 101):
        exact = catalan_number(n-1)
        asymp = 4**n/(4*sp.sqrt(sp.pi)*n**(sp.Rational(3,2)))
        rel_err = abs(exact - asymp)/exact
        if n % 10 == 0:
            print(f"n = {n}: Relative error = {sp.N(rel_err):.5f}")
    # For n=100, the relative error should be less than 1%
    assert sp.N(rel_err) < 0.01, f"Relative error at n=100 is too high: {rel_err}"
    print("Asymptotic estimates test passed.")

# =============================================================================
# Part IV. Generalized Model: Full Ternary Trees and Their Cherry Counts
# =============================================================================

def generate_full_ternary_trees(n):
    """
    Generate all full ternary trees with n leaves.
    In a full ternary tree, every internal node has exactly 3 children.
    Trees exist only for n that can be partitioned into three positive integers.
    """
    if n == 1:
        return ["L"]
    trees = []
    for i in range(1, n-1):
        for j in range(1, n-i):
            k = n - i - j
            if k < 1:
                continue
            left = generate_full_ternary_trees(i)
            middle = generate_full_ternary_trees(j)
            right = generate_full_ternary_trees(k)
            for L in left:
                for M in middle:
                    for R in right:
                        trees.append((L, M, R))
    return trees

def count_cherries_ternary(tree):
    """
    In full ternary trees, define a cherry as an internal node whose all three children are leaves.
    """
    if tree == "L":
        return 0
    if isinstance(tree, tuple) and len(tree) == 3:
        left, mid, right = tree
        is_cherry = 1 if (left == "L" and mid == "L" and right == "L") else 0
        return is_cherry + count_cherries_ternary(left) + count_cherries_ternary(mid) + count_cherries_ternary(right)
    return 0

def test_ternary_trees():
    """Test full ternary trees and cherry counts for n = 1 to 6."""
    print("\nTesting full ternary trees for n = 1 to 6:")
    ternary_data = {}
    for n in range(1, 7):
        trees = generate_full_ternary_trees(n)
        total = len(trees)
        total_cherries = sum(count_cherries_ternary(t) for t in trees)
        freq = {}
        for t in trees:
            c = count_cherries_ternary(t)
            freq[c] = freq.get(c, 0) + 1
        ternary_data[n] = {'trees': total, 'total_cherries': total_cherries, 'frequency': freq}
        print(f"n = {n}: Total trees = {total}, Total cherries = {total_cherries}, Distribution = {freq}")
    return ternary_data

# =============================================================================
# Part V. Unified Test Suite Execution
# =============================================================================

def run_test_suite():
    x, y, z = sp.symbols('x y z')
    max_n = 8  # Use max_n = 8 for binary trees
    
    # 1. Test full binary tree generation.
    print("Testing full binary tree generation:")
    for n in range(1, 10):
        trees = generate_full_binary_trees(n)
        expected = catalan_number(n-1) if n > 1 else 1
        assert len(trees) == expected, f"Tree count mismatch at n={n}: expected {expected}, got {len(trees)}"
    print("Full binary tree generation tests passed.")
    
    # 2. Test invariants for full binary trees.
    print("\nTesting invariants for full binary trees:")
    assert count_cherries("L") == 0, "Leaf should have 0 cherries."
    tree2 = generate_full_binary_trees(2)[0]
    assert count_cherries(tree2) == 1, "n=2 tree should have 1 cherry."
    tree4 = generate_full_binary_trees(4)
    counts = [count_cherries(t) for t in tree4]
    assert counts.count(1) == 4 and counts.count(2) == 1, "n=4 cherry distribution mismatch."
    print("Invariant computation tests passed.")
    
    # 3. Test joint generating function for full binary trees.
    print("\nTesting joint generating function for full binary trees:")
    GF_joint = aggregate_joint_GF(max_n, generate_full_binary_trees)
    print("Joint GF (symbolic form):")
    sp.pprint(GF_joint)
    # Specialize y=1 and z=1 to recover univariate GF.
    GF_univar = sp.simplify(GF_joint.subs({y:1, z:1}))
    Catalan_GF = (1 - sp.sqrt(1-4*x))/2
    # Since our aggregation is only for n=1..max_n, compare the truncated series.
    trunc_GF = sp.series(GF_univar, x, 0, max_n+1).removeO().expand()
    trunc_Cat = sp.series(Catalan_GF, x, 0, max_n+1).removeO().expand()
    diff_GF = sp.simplify(trunc_GF - trunc_Cat)
    print("\nTruncated series expansion difference (GF_univar - Catalan_GF) up to order", max_n, ":")
    sp.pprint(diff_GF)
    assert diff_GF == 0, "Univariate joint GF does not match Catalan GF (truncated series expansion check failed)."
    print("Joint GF tests passed.")
    
    # 4. Test asymptotic estimates for F(x,1).
    test_asymptotic_estimates()
    
    # 5. Test extension: full ternary trees.
    test_ternary_trees()
    
    print("\nAll tests passed successfully.")

if __name__ == "__main__":
    run_test_suite()
