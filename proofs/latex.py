"""A minimal stand-in for `proofs.latex` so these scripts run outside the cbmr-proofs repo.

The real `Proof` renders LaTeX. This one only needs to do the part that makes a proof a proof:
hold a claim's label, check that the expression it rests on is *identically* zero, and fail
loudly if it is not. When these move into cbmr-proofs this file goes away and the import
resolves to the real one; the claim bodies are written against that interface and not this.
"""
import sympy as sp


class ProofFailure(AssertionError):
    """A claim whose expression did not reduce to zero."""


class Proof:
    def __init__(self, name, title, preamble=""):
        self.name = name
        self.title = title
        self.preamble = preamble
        self.claims = []

    def claim(self, label, expression, shown):
        """Record a claim and verify that ``expression`` is identically zero.

        Accepts a scalar or a matrix. A claim that cannot be *shown* zero is a failure rather
        than a warning: these are checked in CI, so a silent pass would be worse than no proof.
        """
        if isinstance(expression, sp.MatrixBase):
            ok = sp.simplify(expression).is_zero_matrix
        elif hasattr(expression, "is_ZeroMatrix") and expression.is_ZeroMatrix:
            ok = True
        else:
            ok = sp.simplify(expression) == 0
        if not ok:
            raise ProofFailure(f"{self.name}: {label!r} did not reduce to zero: {expression}")
        self.claims.append((label, shown))
        return expression

    def report(self):
        print(f"=== {self.title} ({self.name}): {len(self.claims)} claims verified")
        for label, shown in self.claims:
            print(f"  [ok] {label}")
            print(f"       {shown}")
