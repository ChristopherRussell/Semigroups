#############################################################################
##
##  equivlatt.gd
##  Copyright (C) 2019                                    Christopher Russell
##
##  Licensing information can be found in the README file of this package.
##
#############################################################################
##

# Semigroups of uf elements
DeclareCategory("IsUFElement",
                IsAssociativeElement and IsMultiplicativeElementWithInverse and
                IsCommutativeElement and IsIdempotent and IsComponentObjectRep);

BindGlobal("UFElementFamily",
           NewFamily("UFElementFamily",
           IsUFElement));

BindGlobal("UFElementType",
           NewType(UFElementFamily, IsUFElement and IsComponentObjectRep));

# Operations for elements
DeclareOperation("UFFromBlocksNC",
                 [IsList]);
DeclareOperation("UFFromTableNC",
                 [IsList]);
DeclareOperation("UFFromBlocks",
                 [IsList]);
DeclareOperation("UFFromTable",
                 [IsList]);

# Attributes for elements
DeclareAttribute("UFSize",
                 IsUFElement);
DeclareAttribute("UFTable",
                 IsUFElement);
DeclareAttribute("UFBlocks",
                 IsUFElement);
DeclareAttribute("UFNrBlocks",
                 IsUFElement);

# uf semigroups
DeclareCategoryCollections("IsUFElement");
DeclareSynonymAttr("IsUFSemigroup",
                   IsSemigroup and IsUFElementCollection);

#InstallTrueMethod(IsGeneratorsOfEnumerableSemigroup,
#                  IsUFElementCollection);
InstallTrueMethod(IsFinite,
                  IsUFSemigroup);
InstallTrueMethod(IsSemilatticeAsSemigroup,
                  IsUFSemigroup);

# Attributes for uf semigroups
DeclareAttribute("UFSemigroupDegree",
                 IsUFSemigroup);
