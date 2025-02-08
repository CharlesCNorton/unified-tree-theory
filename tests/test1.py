"""
Unified Framework for Tree Invariants in Full Binary Trees

This script performs the following:
1. Generates all full binary trees with n leaves for a range of n.
2. Computes the following invariants for each tree:
   - Cherries: number of internal nodes whose both children are leaves.
   - Colless Index: sum over all internal nodes of |L(v) - R(v)|, where L(v) and R(v) 
     are the numbers of leaves in the left and right subtrees.
   - Sackin Index: sum of the depths of all leaves (with the root at depth 0).
3. Aggregates the results to compute:
   - Total invariant value for each n.
   - Frequency distribution of the invariant values for each n.
4. Constructs univariate generating functions from the brute-force data.
5. Prints detailed tables and (optionally) plots the distributions.
6. This framework is designed to be extendable to other invariants and to support 
   a joint (multivariate) generating function approach in the future.

Authors:
  Charles C. Norton
  OpenAI's o3-mini-high
Date: February 5th, 2025
"""

import sympy as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Part 1: Tree Generation and Invariant Functions
# -----------------------------

def generate_trees(n):
    """
    Recursively generate all full binary trees with n leaves.
    Representation:
      - A leaf is represented by the string "L".
      - An internal node is represented as a tuple (left, right).
    """
    if n == 1:
        return ["L"]
    trees = []
    for i in range(1, n):
        left_subtrees = generate_trees(i)
        right_subtrees = generate_trees(n - i)
        for left in left_subtrees:
            for right in right_subtrees:
                trees.append((left, right))
    return trees

def count_leaves(tree):
    """Return the number of leaves in the tree."""
    if tree == "L":
        return 1
    left, right = tree
    return count_leaves(left) + count_leaves(right)

def count_cherries(tree):
    """
    Count cherries in a full binary tree.
    A cherry is an internal node whose both children are leaves.
    """
    if tree == "L":
        return 0
    left, right = tree
    is_cherry = 1 if (left == "L" and right == "L") else 0
    return is_cherry + count_cherries(left) + count_cherries(right)

def colless_index(tree):
    """
    Compute the Colless index of a full binary tree.
    For each internal node v, with L(v) and R(v) the number of leaves in its left
    and right subtrees respectively, the Colless index is the sum over v of |L(v) - R(v)|.
    """
    if tree == "L":
        return 0
    left, right = tree
    L = count_leaves(left)
    R = count_leaves(right)
    return abs(L - R) + colless_index(left) + colless_index(right)

def sackin_index(tree, depth=0):
    """
    Compute the Sackin index of a full binary tree.
    The Sackin index is the sum of the depths of all leaves (with the root at depth 0).
    """
    if tree == "L":
        return depth
    left, right = tree
    return sackin_index(left, depth + 1) + sackin_index(right, depth + 1)

# -----------------------------
# Part 2: Aggregation Over Tree Invariants
# -----------------------------

def aggregate_invariants(max_n):
    """
    For each n in 1..max_n, generate all full binary trees and compute:
      - The total number of trees (Catalan number).
      - The total cherries, total Colless, and total Sackin values.
      - The frequency distribution of each invariant.
    Returns:
      A dictionary keyed by invariant name, each containing a dictionary keyed by n.
    """
    data = {}
    invariants = ['cherries', 'colless', 'sackin']
    # For each invariant, we record:
    #   - total value (sum over all trees)
    #   - frequency distribution: a dictionary mapping invariant value -> count.
    for inv in invariants:
        data[inv] = {}
    
    for n in range(1, max_n + 1):
        trees = generate_trees(n)
        total_trees = len(trees)
        freq = {inv: {} for inv in invariants}
        total_vals = {inv: 0 for inv in invariants}
        for tree in trees:
            c = count_cherries(tree)
            col = colless_index(tree)
            sack = sackin_index(tree)
            # Update frequency distributions
            for inv, val in zip(invariants, [c, col, sack]):
                freq[inv][val] = freq[inv].get(val, 0) + 1
                total_vals[inv] += val
        # Store results for each invariant
        for inv in invariants:
            data[inv][n] = {'total': total_vals[inv],
                            'frequency': freq[inv],
                            'trees': total_trees}
    return data

# -----------------------------
# Part 3: Construct Univariate Generating Functions from Data
# -----------------------------

def build_generating_function(data_dict, max_n):
    """
    Given a dictionary mapping n to total invariant values (data_dict[n]['total']),
    build the univariate generating function G(x) = sum_{n>=1} (total invariant) x^n.
    Returns a sympy expression for G(x).
    """
    x = sp.symbols('x')
    G = 0
    for n in range(1, max_n + 1):
        G += data_dict[n]['total'] * x**n
    return sp.simplify(G)

# -----------------------------
# Part 4: Main Execution and Results Display
# -----------------------------

max_n = 8  # For brute-force, choose a moderate max_n (Catalan numbers grow rapidly)

print("Generating full binary trees and aggregating invariants for n = 1 to", max_n)
aggregated_data = aggregate_invariants(max_n)

# Print aggregated totals and frequency distributions
for inv in ['cherries', 'colless', 'sackin']:
    print(f"\nInvariant: {inv.capitalize()}")
    for n in range(1, max_n + 1):
        total_val = aggregated_data[inv][n]['total']
        total_trees = aggregated_data[inv][n]['trees']
        freq = aggregated_data[inv][n]['frequency']
        print(f"n = {n} (Total trees = {total_trees}): Total {inv} = {total_val}, Distribution = {freq}")

# Build and print univariate generating functions for each invariant (total value)
G_funcs = {}
for inv in ['cherries', 'colless', 'sackin']:
    G_funcs[inv] = build_generating_function(aggregated_data[inv], max_n)
    print(f"\nUnivariate generating function for total {inv}:")
    sp.pprint(G_funcs[inv])

# -----------------------------
# Part 5: (Optional) Plotting Distributions for a Selected n
# -----------------------------
import matplotlib.pyplot as plt

def plot_distribution(invariant, n, data):
    """Plot the frequency distribution of the specified invariant for trees with n leaves."""
    freq = data[invariant][n]['frequency']
    values = sorted(freq.keys())
    counts = [freq[val] for val in values]
    plt.figure(figsize=(6,4))
    plt.bar(values, counts, color='skyblue')
    plt.xlabel(f"{invariant.capitalize()} Value")
    plt.ylabel("Frequency")
    plt.title(f"Distribution of {invariant.capitalize()} for n = {n} (Total trees = {data[invariant][n]['trees']})")
    plt.xticks(values)
    plt.grid(True, linestyle='--')
    plt.show()

# Plot distributions for n = 7
for inv in ['cherries', 'colless', 'sackin']:
    plot_distribution(inv, 7, aggregated_data)

print("\nUnified framework execution complete. Data has been aggregated and generating functions constructed.")
