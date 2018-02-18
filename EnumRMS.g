###############################################################################
# RMS search - groups and actions
###############################################################################
MatrixTransposeGroup := function()
  trans_group := SymmetricGroup(2);
  trans_act := function(x, g)
    if g = () then
      return x;
    fi;
    return TransposedMat(x);
  end;
  return [trans_group, trans_act];
end;

MatrixOverGroupAutomorphGroup := function(G)
  local auto_grp, auto_act;
  auto_grp := AutomorphismGroup(G);
  auto_act := function(x, g)
    return List(x, row -> List(row, entry -> entry ^ g));
  end;
  return [auto_grp, auto_act];
end;

RMSNormalMatrixColMultGroup := function(G, nr_cols)
  local mult_grp, mult_act;
  mult_grp := SymmetricGroup(nr_cols);
  mult_act := function(x, g);
    return List(x, row -> (row[1 ^ (g ^ -1)] ^ -1) * Permuted(row, g));
  end;
  return [mult_grp, mult_act];
end;

RMSNormalMatrixRowMultGroup := function(G, nr_rows, nr_cols)
  local mult_grp, mult_act;
  mult_grp := SymmetricGroup(nr_rows);
  mult_act := function(x, g)
    local elms;
    elms := List(x[1 ^ (g ^ -1)], a -> a ^ -1);
    return Permuted(List([1 .. nr_rows], i ->
    List([1 .. nr_cols], j -> x[i][j] * elms[j])), g);
  end;
  return [mult_grp, mult_act];
end;

RMSNormalMatrixGroup := function(G, nr_rows, nr_cols)
  local G1, G2, G3, D, p1, p2, p3, act_D;
  G1 := RMSNormalMatrixRowMultGroup(G, nr_rows, nr_cols);
  G2 := RMSNormalMatrixColMultGroup(G, nr_cols);
  G3 := RMSNormalMatrixAutomorphGroup(G, nr_rows, nr_cols);
  D := DirectProduct(G1[1], G2[1], G3[1]);
  p1 := Projection(D, 1);
  p2 := Projection(D, 2);
  p3 := Projection(D, 3);
  act_D := function(x, g)
    return G1[2](G2[2](G3[2](x, g ^ p3), g ^ p2), g ^ p1);
  end;
  return [D, act_D];
end;

RMSNormalMatrixSpace := function(G, rows, cols)
  local col_space, first;
  col_space := Cartesian(Concatenation([Group(())], List([2 .. cols], a -> G)));
  first := Cartesian(List([1 .. cols], a -> Group(())));
  return Cartesian(Concatenation([first],
                                 List([2 .. rows], a -> col_space)));
end;

# RMSNormalMatrixSpace := function(G, nr_rows, nr_cols)
#   local col_space, first;
#   col_space := Cartesian(List([2 .. nr_cols], a -> G));
#   return Cartesian(List([2 .. nr_rows], a -> col_space));
# end;

###############################################################################
# RMS Search
###############################################################################
FindRMSTripleParametersByOrder := function(k)
  local size_G, d, e, out;
  out := [];
  for d in DivisorsInt(k) do
    size_G := d;
    # |I| * |J| = k/d
    for e in DivisorsInt(k / d) do
      # We only want one of (|I|,|J|) = (a,b), (b,a)
      if e <= k / d / e then
        Add(out, [d, e, k / d / e]);
      fi;
    od;
  od;
  return out;
end;

RMSByOrder := function(k)
  local out, parameters, H, p, G;
  out := [];
  parameters := FindRMSTripleParametersByOrder(k);
  for p in parameters do
    for G in AllSmallGroups(p[1]) do
      H := Image(IsomorphismPermGroup(G));
      mats := RMSByParameters(H, p[2], p[3]);
      Append(out, List(mats, mat -> ReesMatrixSemigroup(H, mat)));
    od;
  od;
  return out;
end;

RMSByParameters := function(G, m, n)
  local out, MG, domain, orbs;
  if Size(G) = 1 or m = 1 or n = 1 then
    return [List([1 .. m], i -> List([1 .. n], j -> ()))];
  elif m < n then
    out := RMSByParameters(G, n, m);
    Apply(out, mat -> TransposedMat(mat));
  fi;

  MG := RMSNormalMatrixGroup(G, m, n);
  domain := RMSNormalMatrixSpace(G, m, n);
  orbs := OrbitsDomain(MG[1], domain, MG[2]);
  return List(orbs, orb -> Representative(orb));
end;

# RMSByParameters := function(grp, m, n)
#
#   G := SmallGroup(IdSmallGroup(grp));
#   if Size(G) = 1 or m = 1 or n = 1 then
#     return [List([1 .. m], i -> List([1 .. n], j -> ()))];
#
#   elif m < n then
#     out := RMSByParameters(G, n, m);
#     Apply(out, mat -> TransposedMat(mat));
#
#   elif IsSimpleGroup(G) then # calculate from 1 less rows case
#     MG := RMSNormalMatrixGroup(G, m, n);
#     data := RMSByParameters(G, m - 1, n);
#     extensions := Cartesian([1 .. n], i -> G);
#     candidates := EmptyPlist(Size(data) * Size(extensions));
#     for d in data do
#       Append(candidates, List(extensions, e -> Concatenation(d, [e])));
#     od;
#     return List(Orbits(MG[1], candidates, MG[2]), orb -> Representative(orb));
#
#   else
#     out := [];
#     N := MinimalNormalSubgroup(G);
#     data := RMSByParameters(G/N, m, n);
#     Q := SmallGroup(IdSmallGroup(G/N));
#     mon := IsomorphicSubgroups(G, Q)[1];
#     MG := RMSNormalMatrixGroup(G, m, n);
#     for d in data do
#       candidates := Cartesian(List(d, row -> Cartesian(List(row, g -> N * g ^
#                     mon))));
#       orbs := List(Orbits(MG[1], candidates, MG[2]));
#       Append(out, List(orbs, orb -> Representative(orb)));
#     od;
#     return out;
#   else # decide approach
#     #calculate from maximal homomorphic group image.
#     N := MinimalNormalSubgroup(G);
#
#     #or calculate from 1 less row case
# end;

