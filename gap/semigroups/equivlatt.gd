#############################################################################
##
##  equivlatt.gd
##  Copyright (C) 2019                                    Christopher Russell
##
##  Licensing information can be found in the README file of this package.
##
#############################################################################
##

# Equivalence lattice elements
DeclareCategory("IsEquivalenceLatticeElement",
                IsAssociativeElement and IsMultiplicativeElementWithInverse and
                IsCommutativeElement and IsIdempotent and IsAttributeStoringRep
                and IsPositionalObjectRep);

BindGlobal("EquivalenceLatticeElementFamily",
           NewFamily("EquivalenceLatticeElementFamily",
           IsEquivalenceLatticeElement));

BindGlobal("EquivalenceLatticeElementType",
           NewType(EquivalenceLatticeElementFamily,
           IsEquivalenceLatticeElement and IsPositionalObjectRep));

# Operations for creating equivalence lattice elements
DeclareOperation("EquivalenceLatticeElementNC",
                 [IsList]);
DeclareOperation("EquivalenceLatticeElement",
                 [IsList]);
DeclareOperation("EquivalenceLatticeElement",
                 [IsDigraph]);

# Attributes for equivalence lattice elements
DeclareAttribute("EquivalenceLatticeElementDigraph",
                 IsEquivalenceLatticeElement);
DeclareAttribute("EquivalenceLatticeElementPartition",
                 IsEquivalenceLatticeElement);
DeclareAttribute("EquivalenceLatticeElementDegree",
                 IsEquivalenceLatticeElement);

# Equivalence lattice semigroups
DeclareCategoryCollections("IsEquivalenceLatticeElement");
DeclareSynonymAttr("IsEquivalenceLattice",
                   IsSemigroup and IsEquivalenceLatticeElementCollection);

InstallTrueMethod(IsGeneratorsOfEnumerableSemigroup,
                  IsEquivalenceLatticeElementCollection);
InstallTrueMethod(IsFinite,
                  IsEquivalenceLattice);
InstallTrueMethod(IsSemilatticeAsSemigroup,
                  IsEquivalenceLattice);

# Attributes for equivalence lattice semigroups
DeclareAttribute("EquivalenceLatticeDegree",
                 IsEquivalenceLattice);
