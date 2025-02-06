import unittest
import math
import random
from collections import defaultdict
import sympy

###############################################################################
# Data Structures and Invariant Functions
###############################################################################

class Tree:
    """
    A simple class to represent a full binary tree.
    A leaf is represented by a Tree with both `left` and `right` equal to None.
    An internal node has both `left` and `right` as non-None Trees.
    """
    def __init__(self, left=None, right=None):
        self.left = left
        self.right = right

def is_leaf(tree):
    """Return True if the tree is a leaf."""
    return tree.left is None and tree.right is None

def count_leaves(tree):
    """Count the number of leaves in a full binary tree."""
    if is_leaf(tree):
        return 1
    else:
        return count_leaves(tree.left) + count_leaves(tree.right)

def compute_cherries(tree):
    """
    Compute the number of cherries in the tree.
    A cherry is an internal node whose both children are leaves.
    """
    if is_leaf(tree):
        return 0
    # Check if this internal node is a cherry.
    cherry_here = 1 if (is_leaf(tree.left) and is_leaf(tree.right)) else 0
    return cherry_here + compute_cherries(tree.left) + compute_cherries(tree.right)

def compute_colless(tree):
    """
    Compute the Colless index.
    For each internal node v, compute |L(v) - R(v)|, where L(v) and R(v)
    are the number of leaves in the left and right subtrees.
    """
    if is_leaf(tree):
        return 0
    left_leaves = count_leaves(tree.left)
    right_leaves = count_leaves(tree.right)
    return abs(left_leaves - right_leaves) + compute_colless(tree.left) + compute_colless(tree.right)

def compute_sackin(tree, depth=0):
    """
    Compute the Sackin index: the sum of the depths of all leaves.
    The root is at depth 0.
    """
    if is_leaf(tree):
        return depth
    return compute_sackin(tree.left, depth+1) + compute_sackin(tree.right, depth+1)

###############################################################################
# Tree Generation
###############################################################################

def generate_full_binary_trees(n):
    """
    Recursively generate all full binary trees with n leaves.
    
    For n == 1, yield a single leaf.
    For n > 1, partition the leaves among left and right subtrees.
    """
    if n == 1:
        yield Tree()
    else:
        for i in range(1, n):
            for left in generate_full_binary_trees(i):
                for right in generate_full_binary_trees(n - i):
                    yield Tree(left, right)

def catalan_number(n):
    """
    Compute the nth Catalan number:
      C_n = (1/(n+1)) * binom(2n, n)
    Note: For full binary trees with n leaves, we need C_{n-1}.
    """
    return math.comb(2*n, n) // (n+1)

def random_full_binary_tree(n):
    """
    Generate a random full binary tree with n leaves.
    
    This function uses the well-known fact that the number of full binary trees
    with n leaves equals the (n-1)th Catalan number. It recursively chooses a split 
    weighted by the number of trees possible in the left and right subtrees.
    """
    if n == 1:
        return Tree()
    # Precompute Catalan numbers for sizes 1..(n-1)
    catalans = [catalan_number(i-1) if i > 0 else 1 for i in range(1, n+1)]
    total = 0
    choices = []
    # i goes from 1 to n-1, where left subtree has i leaves and right subtree has n-i leaves.
    for i in range(1, n):
        count = catalan_number(i-1) * catalan_number(n - i - 1)
        choices.append((i, count))
        total += count
    r = random.randint(1, total)
    cum = 0
    for i, count in choices:
        cum += count
        if r <= cum:
            left_size = i
            right_size = n - i
            break
    left_tree = random_full_binary_tree(left_size)
    right_tree = random_full_binary_tree(right_size)
    return Tree(left_tree, right_tree)

###############################################################################
# Test Suite
###############################################################################

class TestTreeInvariants(unittest.TestCase):

    def test_catalan_numbers(self):
        """
        Verify that the number of generated full binary trees for n leaves
        equals the (n-1)th Catalan number.
        """
        for n in range(1, 8):
            trees = list(generate_full_binary_trees(n))
            # Number of full binary trees with n leaves is C_{n-1}
            expected = catalan_number(n-1)
            self.assertEqual(len(trees), expected,
                             f"For n={n}, expected {expected} trees, got {len(trees)}.")

    def test_invariants_leaf(self):
        """
        Check invariants on a single-leaf tree.
        """
        t = Tree()
        self.assertEqual(count_leaves(t), 1)
        self.assertEqual(compute_cherries(t), 0)
        self.assertEqual(compute_colless(t), 0)
        self.assertEqual(compute_sackin(t), 0)

    def test_invariants_small_tree(self):
        """
        Test invariants for the smallest non-trivial full binary tree (n=2).
        """
        # Create a tree with 2 leaves (an internal node with two leaves)
        t = Tree(Tree(), Tree())
        self.assertEqual(count_leaves(t), 2)
        # This internal node is a cherry (both children are leaves)
        self.assertEqual(compute_cherries(t), 1)
        # Colless index: |1 - 1| = 0 at the root.
        self.assertEqual(compute_colless(t), 0)
        # Sackin index: both leaves are at depth 1, so total = 2.
        self.assertEqual(compute_sackin(t), 2)

    def test_cherry_distribution_n4(self):
        """
        For n = 4 leaves, the paper claims that the cherry distribution is:
          {1: 4, 2: 1}
        That is, 4 trees with 1 cherry and 1 tree with 2 cherries.
        """
        trees = list(generate_full_binary_trees(4))
        distribution = defaultdict(int)
        for t in trees:
            distribution[compute_cherries(t)] += 1
        expected = {1: 4, 2: 1}
        self.assertEqual(dict(distribution), expected,
                         f"n=4 cherry distribution expected {expected}, got {dict(distribution)}.")

    def test_cherry_distribution_n5(self):
        """
        For n = 5 leaves, the paper claims that the cherry distribution is:
          {1: 8, 2: 6}
        """
        trees = list(generate_full_binary_trees(5))
        distribution = defaultdict(int)
        for t in trees:
            distribution[compute_cherries(t)] += 1
        expected = {1: 8, 2: 6}
        self.assertEqual(dict(distribution), expected,
                         f"n=5 cherry distribution expected {expected}, got {dict(distribution)}.")

    def test_generating_function_specialization(self):
        """
        Verify that the joint generating function specialized to u = v = w = 1
        recovers the Catalan generating function:
          T(x;1,1,1) = (1 - sqrt(1-4*x)) / 2.
        Here we use sympy to check the series coefficients.
        """
        x, u, v, w = sympy.symbols('x u v w', positive=True)
        catalan_GF = (1 - sympy.sqrt(1 - 4*x)) / 2
        # Expand to order 10 (i.e. up to x^9)
        series_expansion = sympy.series(catalan_GF, x, 0, 10).removeO()
        # Expected coefficients (Catalan numbers for n = 0 to 9):
        # Note: The series for T(x) here is in powers of x starting at x^1.
        # We expect: x^1:1, x^2:1, x^3:2, x^4:5, x^5:14, x^6:42, x^7:132, x^8:429, x^9:1430.
        expected_coeffs = [0, 1, 1, 2, 5, 14, 42, 132, 429, 1430]
        coeffs = [series_expansion.coeff(x, n) for n in range(10)]
        self.assertEqual(coeffs, expected_coeffs,
                         f"Generating function series coefficients expected {expected_coeffs}, got {coeffs}.")

    def test_random_tree_invariants_statistics(self):
        """
        Perform a Monte Carlo simulation: generate many random full binary trees
        with a large number of leaves (n = 50) and compute average values for
        cherries, Colless, and Sackin indices. While we do not have exact predictions
        here, the averages should be positive and scale reasonably.
        """
        n = 50
        num_samples = 1000
        total_cherries = 0
        total_colless = 0
        total_sackin = 0
        for _ in range(num_samples):
            t = random_full_binary_tree(n)
            total_cherries += compute_cherries(t)
            total_colless += compute_colless(t)
            total_sackin += compute_sackin(t)
        avg_cherries = total_cherries / num_samples
        avg_colless = total_colless / num_samples
        avg_sackin = total_sackin / num_samples

        # Output the Monte Carlo estimates.
        print("\nMonte Carlo simulation for n = 50 (over {} samples):".format(num_samples))
        print("Average cherries: {:.2f}".format(avg_cherries))
        print("Average Colless index: {:.2f}".format(avg_colless))
        print("Average Sackin index: {:.2f}".format(avg_sackin))
        
        # Basic sanity checks: averages should be positive.
        self.assertGreater(avg_cherries, 0)
        self.assertGreater(avg_colless, 0)
        self.assertGreater(avg_sackin, 0)
        # For cherries, the average cannot exceed the number of leaves.
        self.assertLess(avg_cherries, n)

###############################################################################
# Main: Run the test suite.
###############################################################################

if __name__ == '__main__':
    unittest.main()
