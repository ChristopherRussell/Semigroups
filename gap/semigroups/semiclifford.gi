InstallMethod(SemilatticeDigraphToSemilatticeSemigroup,
"for a digraph",
[IsDigraph],
function(digraph)
  local red, top, im, max, n, S, x, i, j;

  if IsMeetSemilatticeDigraph(digraph) then
    return SemilatticeDigraphToSemilatticeSemigroup(Reversed(digraph))
  elif not IsJoinSemilatticeDigraph(digraph) then
    ErrorNoReturn("Semigroups: SemilatticeDigraphToSemilatticeSemigroup ",
                  "usage,\n", "the argument must be a join semilattice ",
                  "digraph or a meet semilattice digraph,");
  fi;

  red := DigraphReflexiveTransitiveReduction(digraph);
  top := DigraphTopologicalSort(digraph);
  # im[i] will store the image of the idempotent partial perm corresponding to
  # vertex i of the arugment <digraph>
  im         := [];
  im[top[1]] := [];
  max        := 1;

  n := DigraphNrVertices(digraph);
  # For each vertex, the corresponding idempotent has an image
  # containing all images of idempotents below it.
  for i in [2 .. n] do
    im[top[i]] := Union(List(OutNeighboursOfVertex(red, top[i]), j -> im[j]));
    # When there is only one neighbour, we must add a point to the image to
    # distinguish the two idempotent partial perms.
    if Length(OutNeighboursOfVertex(red, top[i])) = 1 then
      Add(im[top[i]], max);
      max := max + 1;
    fi;
  od;

  # Add the generators one by one
  S := Semigroup(PartialPerm(im[top[n]], im[top[n]]));
  for j in Reversed(top){[2 .. n]} do
    x := PartialPerm(im[j], im[j]);
    if x in S then
      continue;
    fi;
    S := ClosureSemigroup(S, x);
    if Size(S) = n then
      return S;
    fi;
  od;

  return S;
end);

InstallMethod(CliffordSemigroupConstructor,
"for a digraph, list and list",
[IsDigraph, IsList, IsList],
function(digraph, groups, homs)
  local n, red, nr_edge, edges, top, im, max, S, x, i, j;

  if not IsJoinSemilatticeDigraph(digraph) then
    ErrorNoReturn("the first argument should be a join semilattice");
  fi;

  n       := DigraphNrVertices(digraph);
  red     := DigraphReflexiveTransitiveReduction(digraph);
  nr_edge := DigraphNrEdges(red);
  edges   := DigraphEdges(red);
  if Length(groups) <> n then
    ErrorNoReturn("there should be a group corresponding to each vertex");
  elif not ForAll(groups, IsPermGroup) then
    ErrorNoReturn("the second argument should be a list of perm groups");
  elif Length(homs) <> nr_edge(red) then
    ErrorNoReturn("there should be a homomorphism corresponding to each edge");
  elif not ForAll(homs, IsGroupHomomorphism) then
    ErrorNoReturn("the third argument should be a list of group homomorphism");
  elif not ForAll([1 .. nr_edge], i -> Source(hom[i]) = groups(edges[i][1]) and
                                  Range(hom[i]) = groups(edges[i][2])) then
    ErrorNoReturn("homomorphisms must match the vertex groups joined by edges");
  fi;

  top          := DigraphTopologicalSort(red);
  deg          := NrMovedPoints(groups[top[1]]);
  # preim_doms are the disjoint domains of the groups, which have size equal to
  # the degree of the group
  preim_doms         := [];
  preim_doms[top[1]] := [1 .. NrMovedPoints(groups[top[1]])];
  # im_doms are the domains of images of the groups in the partial perm
  # representation, which must be larger than all groups below them in the
  # semilattice
  im_doms      := ShallowCopy(preim_doms);
  for i in [2 .. n] do
    preim_doms[top[i]] := [1 .. NrMovedPoints(groups[top[i]])]
                     + Maximum(preim_doms[top[i - 1]]);
    im_doms[top[i]] := Union(List(OutNeighboursOfVertex(red, top[i]),
                             j -> im_doms[j]));
    Append(im_doms[top[i]], preim_doms[top[i]]);
  od;

  gens := [];
  Append(gens, Generators(groups[top[1]]));
  for i in [2 .. n] do
    u   := top[i];
    G   := groups[u];

    for g in Generators(G) do
      im := ShallowCopy(im_doms[u]);
      for v in OutNeighboursOfVertex(red, u) do
        x := g ^ homs[Position(edges, [u, v])];
        min := Minimum(preim_doms) - 1;
        for a in im do
          if a in preim_doms[u] then
            a := (a - min) ^ x + min;
          fi;
        od;
      od;
      Add(gens, PartialPerm(im_doms[u], s));
    od;
  od;

  return S;
end);
