TODO

* remove C++11
* gmp problems?
* use namespace in the C++ code
* clean up the memory usage in semigroups.h, interface.cc
* add _report for all methods that call enumerate in semigroups.h
* DClass(TransformationSemigroup, BooleanMat) sometimes seg faults or returns an answers which it shouldn't!
* reporting in semigroups.h (toggle_report and _report)

Gaplint:

* NrLeftBlocks, NrBlocks are in gaplint_ignore since the every way to create 
  a bipartition sets these attribute just after creation. Hence if I install a
  method for these, then it will never be run, hence hurting our code coverage. 

Code coverage:

* bipartition.tst tests factor.gi to 80%

General notes:

* need to check usage of IsSemigroup and IsFinite, the max-plus matrix
  semigroup examples can finite or infinite, we should check if they are finite
  (somehow) and then use RightCayleyGraphSemigroup...

* HasGeneratorsOfMagmaIdeal has been changed to HasGeneratorsOfSemigroupIdeal
  in lots of places.


New features in this branch:

* IsSemigroupWithInverseOp 
* SemigroupOptions is an attribute, has a component generic which is set to
  true or false