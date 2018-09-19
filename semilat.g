SemilatticeDigraphToSemilatticeSemigroup := function(digraph)
  local red, top, im, max, n, S, x, i, j;

  if not IsJoinSemilatticeDigraph(digraph) then
    return fail;
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
end;
