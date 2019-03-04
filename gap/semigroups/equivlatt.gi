# Semigroups
InstallMethod(UFSemigroupDegree, "for a uf semigroup",
[IsUFSemigroup],
function(S)
  return UF_SIZE(S.1!.uf);
end);

InstallMethod(OneImmutable, "for a uf semigroup",
[IsUFSemigroup],
function(S)
  return Objectify(UFElementType, rec(uf := UF_NEW(UFSemigroupDegree(S))));
end);

InstallMethod(String, "for a uf semigroup",
[IsUFSemigroup],
function(S)
  return Concatenation("Semigroup(", String(GeneratorsOfSemigroup(S)));
end);

InstallMethod(PrintObj, "for a uf semigroup",
[IsUFSemigroup],
function(S)
  Print(String(S));
  return;
end);

InstallMethod(ViewString, "for a uf semigroup",
[IsUFSemigroup],
function(S)
  return Concatenation("<uf semigroup of degree ",
                       String(UFSemigroupDegree(S)), ">");
end);

# Elements

InstallMethod(UFFromBlocksNC,
"for a list",
[IsList],
function(b)
  local n, u, block, x;
  n := Sum(b, Length);
  u := UF_NEW(n);
  for block in b do
    for x in block{[2 .. Length(block)]} do
      UF_UNION(u, [block[1], x]);
    od;
  od;
  return Objectify(UFElementType, rec(uf := u, blocks := b));
end);

InstallMethod(UFFromBlocks,
"for a list",
[IsList],
function(b)
  if not ForAll(b, block -> IsList(block) and ForAll(block, IsInt)) then
    ErrorNoReturn("Semigroups: EquivalenceLatticeElement: usage,\n",
                  "the argument should be a list of ints,");
  fi;
  return UFFromBlocksNC(b);
end);

InstallMethod(UFFromTableNC,
"for a list",
[IsList],
function(t)
  local n, u, x;
  n := Length(t);
  u := UF_NEW(n);
  for x in [1 .. n] do
    UF_UNION(u, [x, t[x]]);
  od;
  return Objectify(UFElementType, rec(uf := u, table := t));
end);

InstallMethod(UFFromTable,
"for a list",
[IsList],
function(t)
  if not ForAll(t, IsInt) then
    ErrorNoReturn("Semigroups: EquivalenceLatticeElement: usage,\n",
                  "the argument should be a list of ints,");
  fi;
  return UFFromTableNC(t);
end);

InstallMethod(UFSize,
"for a list",
[IsUFElement],
function(uf)
  return UF_SIZE(uf!.uf);
end);

InstallMethod(UFTable,
"for a list",
[IsUFElement],
function(uf)
  if IsBound(uf!.table) then
    return uf!.table;
  fi;
  return UF_TABLE(uf!.uf);
end);

InstallMethod(UFBlocks,
"for a list",
[IsUFElement],
function(uf)
  return UF_BLOCKS(uf!.uf);
end);

InstallMethod(UFNrBlocks,
"for a list",
[IsUFElement],
function(uf)
  return UF_NR_BLOCKS(uf!.uf);
end);

InstallMethod(String, "for a uf object",
[IsUFElement],
function(uf)
  return String(UFBlocks(uf));
end);

InstallMethod(ViewString, "for a uf object",
[IsUFElement],
function(uf)
  return ViewString(UFBlocks(uf));
end);

InstallMethod(\=, "for two uf objects",
IsIdenticalObj,
[IsUFElement, IsUFElement],
function(uf1, uf2)
  if not UFSize(uf1) = UFSize(uf2) then
    ErrorNoReturn("Semigroups: \=, usage:",
                  "the two uf elements must have the same size,");
  fi;
  return ForAll([1 .. UFSize(uf1)],
                i -> UF_FIND(uf1!.uf, i) = UF_FIND(uf2!.uf, i));
end);

InstallMethod(\*, "for two uf objects",
[IsUFElement, IsUFElement],
function(uf1, uf2)
  if not UFSize(uf1) = UFSize(uf2) then
    ErrorNoReturn("Semigroups: \*, usage:",
                  "the two uf elements must have the same size,");
  fi;
  return Objectify(UFElementType, rec(uf := UF_JOIN(uf1!.uf, uf2!.uf)));
end);

InstallMethod(\<, "for two uf objects",
[IsUFElement, IsUFElement],
function(uf1, uf2)
  local i;
  if not UFSize(uf1) = UFSize(uf2) then
    ErrorNoReturn("Semigroups: \<, usage:",
                  "the two uf elements must have the same size,");
  fi;
  for i in [1 .. UFSize(uf1)] do
    if UF_FIND(uf1!.uf, i) < UF_FIND(uf2!.uf, i) then
      return true;
    elif UF_FIND(uf1!.uf, i) > UF_FIND(uf2!.uf, i) then
      return false;
    fi;
  od;
  return false;
end);

InstallMethod(InverseOp, "for a uf object",
[IsUFElement],
function(uf)
  return uf;
end);

InstallMethod(\^, "for a uf object and a negative int",
[IsUFElement, IsNegInt],
function(uf, i)
  return uf;
end);

InstallMethod(LeftOne, "for a uf object",
[IsUFElement],
function(uf)
  return uf;
end);

InstallMethod(RightOne, "for a uf object",
[IsUFElement],
function(uf)
  return uf;
end);

InstallMethod(RightOne, "for a uf object",
[IsUFElement],
function(uf)
  return uf;
end);

InstallMethod(IsGreensDGreaterThanFunc, "for a uf semigroup",
[IsUFSemigroup],
function(S)
  return {x, y} -> x * y = y;
end);
