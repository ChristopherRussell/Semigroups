# Semigroups
InstallMethod(EquivalenceLatticeDegree, "for a equivalence lattice",
[IsEquivalenceLattice],
function(S)
  return EquivalenceLatticeElementDegree(S.1);
end);

InstallMethod(OneImmutable, "for a equivalence lattice",
[IsEquivalenceLattice],
function(S)
  return EquivalenceLatticeElement(
    List([1 .. EquivalenceLatticeDegree(S)], i -> [i]));
end);

InstallMethod(String, "for a equivalence lattice",
[IsEquivalenceLattice],
function(S)
  return Concatenation("Semigroup(", String(GeneratorsOfSemigroup(S)));
end);

InstallMethod(PrintObj, "for a equivalence lattice",
[IsEquivalenceLattice],
function(S)
  Print(String(S));
  return;
end);

InstallMethod(ViewString, "for a equivalence lattice",
[IsEquivalenceLattice],
function(S)
  return Concatenation("<equivalence lattice of degree ",
                       String(EquivalenceLatticeDegree(S)), ">");
end);

# Elements

InstallMethod(EquivalenceLatticeElement,
"for a list",
[IsList],
function(p)
  if not ForAll(p, IsInt) then
    ErrorNoReturn("Semigroups: EquivalenceLatticeElement: usage,\n",
                  "the argument should be a list of ints,");
  fi;
  return Objectify(EquivalenceLatticeElementType, p);
end);

InstallMethod(EquivalenceLatticeElementNC,
"for a list",
[IsList],
function(p)
  return Objectify(EquivalenceLatticeElementType, p);
end);

InstallMethod(EquivalenceLatticeElement,
"for a digraph",
[IsDigraph],
function(digraph)
  # if not IsReflexiveDigraph(digraph) and IsTransitiveDigraph(digraph) and
  #   IsSymmetricDigraph(digraph) then
  #   ErrorNoReturn("Semigroups: EquivalenceLatticeElement: usage,\n",
  #                 "the argument should be a reflexive, transitive and ",
  #                 "symmetric digraph,");
  # fi;
  return Objectify(EquivalenceLatticeElementType, rec(digraph := digraph));
end);

InstallMethod(EquivalenceLatticeElementDigraph,
"for a equivalence lattice element",
[IsEquivalenceLatticeElement],
function(x)
  local edges, part, i;
  if IsBound(x!.digraph) then
    return x!.digraph;
  fi;
  edges := [];
  for part in EquivalenceLatticeElementPartition(x) do
    for i in [1 .. Length(part) - 1] do
        AddSet(edges, [part[i], part[i + 1]]);
    od;
    AddSet(edges, [part[Length(part)], part[1]]);
  od;
  x!.digraph := DigraphByEdges(edges);
  return x!.digraph;
end);

InstallMethod(EquivalenceLatticeElementPartition,
"for a equivalence lattice element",
[IsEquivalenceLatticeElement],
function(x)
  if IsBound(x) then
    return x;
  fi;
  x!.partition := DigraphStronglyConnectedComponents(
    EquivalenceLatticeElementDigraph(x))!.comps;
  return x!.partition;
end);

InstallMethod(EquivalenceLatticeElementDegree,
"for a equivalence lattice element",
[IsEquivalenceLatticeElement],
function(x)
  if IsBound(x!.partition) then
    return Size(Union(EquivalenceLatticeElementPartition(x)));
  else
    return DigraphNrVertices(EquivalenceLatticeElementDigraph(x));
  fi;
end);

InstallMethod(String, "for a equivalence lattice element rep",
[IsEquivalenceLatticeElement],
function(x)
  return String(EquivalenceLatticeElementPartition(x));
end);

InstallMethod(ViewString, "for a equivalence lattice element rep",
[IsEquivalenceLatticeElement],
function(x)
  return ViewString(EquivalenceLatticeElementPartition(x));
end);

InstallMethod(\=, "for two equivalence lattice element reps",
IsIdenticalObj,
[IsEquivalenceLatticeElement, IsEquivalenceLatticeElement],
function(x, y)
  return EquivalenceLatticeElementPartition(x)
         = EquivalenceLatticeElementPartition(y);
end);

InstallMethod(\*, "for two equivalence lattice element reps",
[IsEquivalenceLatticeElement, IsEquivalenceLatticeElement],
function(x, y)
  local Ex, Ey, di, xy;
  if not EquivalenceLatticeElementDegree(x) =
         EquivalenceLatticeElementDegree(y) then
    ErrorNoReturn("Semigroups: \*, usage:",
                  "the two elements must have the same degree,");
  fi;
  Ex := DigraphEdges(EquivalenceLatticeElementDigraph(x));
  Ey := DigraphEdges(EquivalenceLatticeElementDigraph(y));
 # di := DigraphReflexiveTransitiveClosure(DigraphByEdges(Union(Ex, Ey)));
 # return EquivalenceLatticeElement(DigraphStronglyConnectedComponents(di)!.comps);
  
  di := DigraphByEdges(Union(Ex, Ey));
  xy := EquivalenceLatticeElementNC(
    DigraphStronglyConnectedComponents(di).comps);
  SetEquivalenceLatticeElementDegree(xy, EquivalenceLatticeElementDegree(x));
  return xy;
end);

InstallMethod(\<, "for two equivalence lattice element reps",
[IsEquivalenceLatticeElement, IsEquivalenceLatticeElement],
function(x, y)
  if not EquivalenceLatticeElementDegree(x) =
         EquivalenceLatticeElementDegree(y) then
    ErrorNoReturn("Semigroups: \*, usage:",
                  "the two elements must have the same degree,");
  fi;
  return ForAll(EquivalenceLatticeElementPartition(y),
                p -> ForAny(EquivalenceLatticeElementPartition(x),
                            q -> IsSubset(q, p)));
end);

InstallMethod(InverseOp, "for a equivalence lattice element rep",
[IsEquivalenceLatticeElement],
function(x)
  return x;
end);

InstallMethod(\^, "for a equivalence lattice element and a negative int",
              [IsEquivalenceLatticeElement, IsNegInt],
function(x, i)
  return x;
end);

InstallMethod(LeftOne, "for a equivalence lattice element rep",
[IsEquivalenceLatticeElement],
function(x)
  return x;
end);

InstallMethod(RightOne, "for a equivalence lattice element rep",
[IsEquivalenceLatticeElement],
function(x)
  return x;
end);
