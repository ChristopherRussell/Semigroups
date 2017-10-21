DeclareAttribute("MaximumIdempotentSeparatingCongruence", IsInverseSemigroup);
DeclareAttribute("LatticeOfIdempotentSeparatingCongruences",
                 IsInverseSemigroup);
DeclareAttribute("IdempotentSeparatingMonomorphismToMunnSemigroup",
                  IsInverseSemigroup);
DeclareAttribute("FundamentalInverseSemgroupExtensionRepresentation",
  IsInverseSemigroup);

DeclareProperty("IsFundamentalInverseSemigroup", IsSemigroup);
DeclareProperty("IsIdempotentSeparatingCongruence", IsSemigroupCongruence);

DeclareOperation("InverseSemigroupIdempotentSeparatingCongruence",
  [IsInverseSemigroup, IsInverseSemigroup]);
DeclareOperation("InverseSemigroupIdempotentSeparatingCongruenceByKernel",
  [IsInverseSemigroup]);
DeclareOperation(
  "InverseSemigroupIdempotentSeparatingCongruenceByGeneratingPair",
  [IsInverseSemigroup]);

DeclareOperation("InverseSemigroupOrderIdeal", [IsInverseSemigroup, IsList]);
DeclareOperation("InverseSemigroupOrderIdeal",
                  [IsInverseSemigroup, IsMultiplicativeElementWithInverse]);
