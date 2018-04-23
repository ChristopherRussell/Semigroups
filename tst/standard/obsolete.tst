#############################################################################
##
#W  standard/obsolete.tst
#Y  Copyright (C) 2016                                   James D. Mitchell
##
##  Licensing information can be found in the README file of this package.
##
#############################################################################
##
gap> START_TEST("Semigroups package: standard/obsolete.tst");
gap> LoadPackage("semigroups", false);;

#
gap> SEMIGROUPS.StartTest();

# RandomTransformationSemigroup
gap> RandomTransformationSemigroup(2, 2);;
#I  `RandomTransformationSemigroup` is no longer supported
#I  use `RandomSemigroup(IsTransformationSemigroup, 2, 2)` instead!

# RandomTransformationMonoid
gap> RandomTransformationMonoid(2, 2);;
#I  `RandomTransformationMonoid` is no longer supported
#I  use `RandomMonoid(IsTransformationMonoid, 2, 2)` instead!

# RandomPartialPermSemigroup
gap> RandomPartialPermSemigroup(2, 2);;
#I  `RandomPartialPermSemigroup` is no longer supported
#I  use `RandomSemigroup(IsPartialPermSemigroup, 2, 2)` instead!

# RandomPartialPermMonoid
gap> RandomPartialPermMonoid(2, 2);;
#I  `RandomPartialPermMonoid` is no longer supported
#I  use `RandomMonoid(IsPartialPermMonoid, 2, 2)` instead!

# RandomMatrixSemigroup
gap> RandomMatrixSemigroup(2, 2);;
#I  `RandomMatrixSemigroup` is no longer supported
#I  use `RandomSemigroup(IsMatrixOverFiniteFieldSemigroup, 2, 2)` instead!

# RandomMatrixMonoid
gap> RandomMatrixMonoid(2, 2);;
#I  `RandomMatrixMonoid` is no longer supported
#I  use `RandomMonoid(IsMatrixOverFiniteFieldMonoid, 2, 2)` instead!

# DotDClasses 
gap> S := RandomPartialPermMonoid(2, 2);;
#I  `RandomPartialPermMonoid` is no longer supported
#I  use `RandomMonoid(IsPartialPermMonoid, 2, 2)` instead!
gap> DotDClasses(S);;
#I  `DotDClasses` is no longer supported
#I  use `DotString` instead!
gap> DotDClasses(S, rec());;
#I  `DotDClasses` is no longer supported
#I  use `DotString` instead!

# PartialTransformationSemigroup
gap> PartialTransformationSemigroup(5);;
#I  `PartialTransformationSemigroup` is no longer supported
#I  use `PartialTransformationMonoid(5)` instead!

# AsPartialPermSemigroup 
gap> AsPartialPermSemigroup(DualSymmetricInverseMonoid(3));;
#I  `AsPartialPermSemigroup` is no longer supported
#I  use `AsSemigroup(IsPartialPermSemigroup, S)` instead!

# AsPartialPermMonoid 
#gap> AsPartialPermMonoid(FullTransformationMonoid(5));;
# AsBipartitionSemigroup
gap> AsBipartitionSemigroup(FullTransformationMonoid(5));;
#I  `AsBipartitionSemigroup` is no longer supported
#I  use `AsSemigroup(IsBipartitionSemigroup, S)` instead!

# AsBlockBijectionSemigroup
gap> AsBlockBijectionSemigroup(DualSymmetricInverseMonoid(3));;
#I  `AsBlockBijectionSemigroup` is no longer supported
#I  use `AsSemigroup(IsBlockBijectionSemigroup, S)` instead!

# AsTransformationSemigroup
gap> AsTransformationSemigroup(DualSymmetricInverseMonoid(3));;
#I  `AsTransformationSemigroup` is no longer supported
#I  use `AsSemigroup(IsTransformationSemigroup, S)` instead!

# AsMatrixSemigroup
gap> AsMatrixSemigroup(DualSymmetricInverseMonoid(3));;
#I  `AsMatrixSemigroup` is no longer supported
#I  use `AsSemigroup(IsMatrixOverFiniteFieldSemigroup, S)` instead!

# IsomorphismBipartitionSemigroup
gap> IsomorphismBipartitionSemigroup(SymmetricInverseMonoid(3));;
#I  `IsomorphismBipartitionSemigroup` is no longer supported
#I  use `IsomorphismSemigroup(IsBipartitionSemigroup, S)` instead!

# IsomorphismBlockBijectionSemigroup
gap> IsomorphismBlockBijectionSemigroup(SymmetricInverseMonoid(3));;
#I  `IsomorphismBlockBijectionSemigroup` is no longer supported
#I  use `IsomorphismSemigroup(IsBlockBijectionSemigroup, S)` instead!

# IsomorphismMatrixSemigroup
gap> IsomorphismMatrixSemigroup(SymmetricInverseMonoid(3));;
#I  `IsomorphismMatrixSemigroup` is no longer supported
#I  use `IsomorphismSemigroup(IsMatrixOverFiniteFieldSemigroup, S)` instead!

# IsSynchronizingTransformationCollection
gap> IsSynchronizingTransformationCollection([
>  Transformation([2, 3, 4, 1]),
>  Transformation([2, 1]),
>  Transformation([1, 2, 3, 1])], 4);;
#I  `IsSynchronizingTransformationCollection` is no longer supported
#I  use `IsSynchronizingSemigroup(Semigroup(coll), n)` instead!

# FactorisableDualSymmetricInverseSemigroup
gap> FactorisableDualSymmetricInverseSemigroup(3);;
#I  `FactorisableDualSymmetricInverseSemigroup` is no longer supported
#I  use `FactorisableDualSymmetricInverseMonoid(3)` instead!

# SingularFactorisableDualSymmetricInverseSemigroup
gap> SingularFactorisableDualSymmetricInverseSemigroup(3);;
#I  `SingularFactorisableDualSymmetricInverseSemigroup` is no longer supported
#I  use `SingularFactorisableDualSymmetricInverseMonoid(3)` instead!

#
gap> SEMIGROUPS.StopTest();
gap> STOP_TEST("Semigroups package: standard/obsolete.tst");
