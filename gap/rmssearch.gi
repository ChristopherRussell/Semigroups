# We find the possible sizes |G|, |I|, |J| such that |RMS(G, I, J, P)|. We call
# this the 'shape' of (G, I, J)
FindRMSTripleShapeByOrder := function(n)
  local size_G, d, e, out;
  out := [];
  for d in DivisorsInt(n) do
    size_G := d;
    # |I| * |J| = n/d
    for e in DivisorsInt(n/d) do
      # We only want one of (|I|,|J|) = (a,b), (b,a)
      if e <= n/d/e then
        Add(out, [d, e, n/d/e]);
      fi;
    od;
  od;
  return out;
end;

# This method takes as argument a collection of RMS reduces it to a list which
# is unique up to isomorphism.
IsomorphismTestRMSCollection := function(coll)
  local i, j, m, n, copy, out;

  copy := StructuralCopy(coll);
  out := [];
  for i in [1 .. Length(copy)] do
    if IsBound(copy[i]) then
      Add(out, copy[i]);
    fi;
    for j in [Length(copy), Length(copy) - 1 .. i + 1] do
      if IsIsomorphicSemigroup(copy[i], copy[j]) then
        Remove(copy, j);
      fi;
    od;
  od;
  return out;
end;

# This method takes as argument a collection of RMS returns a lists of lists
# which are equivalence classes up to isomorphism.
IsomorphismTestRMSCollectionExperiment := function(coll)
  local i, j, m, n, copy, count, out;

  copy := StructuralCopy(coll);
  out := [];
  count := 1;
  for i in [1 .. Length(copy)] do
    if IsBound(copy[i]) then
      Add(out, [copy[i]]);
    fi;
    for j in [Length(copy), Length(copy) - 1 .. i + 1] do
      if IsIsomorphicSemigroup(copy[i], copy[j]) then
        Add(out[count], copy[j]);
        Remove(copy, j);
      fi;
    od;
    count := count + 1;
  od;
  return out;
end;

###############################################################################
# Finding RMS by calculating orbits over normal matrix space
###############################################################################

RMSNormalMatrixAutomorphGroup := function(G, rows, cols)
  local auto_act, auto_grp;
  auto_grp := AutomorphismGroup(G);
  auto_act := function(x ,g)
    return List(x, row -> List(row, entry -> entry ^ g));
  end;
  return [auto_grp, auto_act];
end;

RMSNormalMatrixColMultGroup := function(G, cols)
  local mult_grp, mult_act;
  mult_grp := SymmetricGroup(cols);
  mult_act := function(x, g);
    return List(x, row -> Permuted(row, g) * row[1 ^ (g ^ -1)] ^ -1);
  end;
  return [mult_grp, mult_act];
end;

RMSNormalMatrixRowMultGroup := function(G, rows, cols)
  local mult_grp, mult_act;
  mult_grp := SymmetricGroup(rows);
  mult_act := function(x, g)
    local elms;
    elms := List(x[1 ^ (g ^ -1)], a -> a ^ -1); 
    return Permuted(List([1 .. rows], i ->
    List([1 .. cols], j -> x[i][j] * elms[j])), g);
  end;
  return [mult_grp, mult_act];
end;

RMSNormalMatrixGroup := function(G, rows, cols)
  local G1, G2, G3, D, p1, p2, p3, act_D;
  G1 := RMSNormalMatrixRowMultGroup(G, rows, cols);
  G2 := RMSNormalMatrixColMultGroup(G, cols);
  G3 := RMSNormalMatrixAutomorphGroup(G, rows, cols);
  D := DirectProduct(G1[1], G2[1], G3[1]);
  p1 := Projection(D, 1);
  p2 := Projection(D, 2);
  p3 := Projection(D, 3);
  # iso := IsomorphismPermGroup(G3[1]);
  # G3 := Image(iso);
  # p3 := CompositionMapping(p3, InverseGeneralMapping(iso));
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

FindRMSByOrderImproved := function(n)
  local s, G, H, mat_grp, mats, orbs, out, shapes;

  shapes := FindRMSTripleShapeByOrder(n);
  out := [];
  for s in shapes do
    Print("\c Current order = ", n, "; and current shape = ", s, "\r");
    for G in AllSmallGroups(s[1]) do
      H := Image(IsomorphismPermGroup(G));
      if s[2] = 1 then
        Add(out, [ReesMatrixSemigroup(H, [List([1 .. s[3]], a -> ())])]);
      elif s[3] = 1 then
        Add(out, [ReesMatrixSemigroup(H, [List([1 .. s[2]], a -> [()])])]);
      elif s[1] = 1 then
        Add(out, [ReesMatrixSemigroup(H, List([1 .. s[2]], a -> List([1 ..
          s[3]], b -> ())))]);
      else
        mats := RMSNormalMatrixSpace(H, s[2], s[3]);
        mat_grp := RMSNormalMatrixGroup(H, s[2], s[3]);
        orbs := OrbitsDomain(mat_grp[1], mats, mat_grp[2]);
        Add(out, List(orbs, orb -> ReesMatrixSemigroup(H, orb[1])));
      fi;
    od;
  od;
  return out;
end;

# Numbers of RMS by order, up to isomorphism except in square matrix case.
# [ 1, 2, 2, 5, 2, 6, 2, 11, 5, 6, 2, 16, 2, 6, 5, 28, 2, 16, 2, 16, 6, 6, 2,
# 42, 5, 6, 11, 15, 2, 19, 2, 85, 5, 6, 5, 47, 2, 6, 6, 41, 2, 22, 2, 15, 12,
# 6, 2, 123, 5, 16, 5, 16, 2, 42, 6, 39, 6, 6, 2, 57, 2, 6, 15 ]

###############################################################################
# RZMS Search
###############################################################################

RZMSMatrixShapesByParametersQuicker := function(m, n)
  local bt, first, next, output, i, act, row_space, size_row_space, G, out,
  shapes, orbs;

  if m = 1 and n = 1 then
    return Set([[0]], [[1]]);
  elif m = 1 and n > 1 then
    return Set([0 .. n], a -> [Concatenation(List([1 .. a], b -> 0),
                                              List([a + 1 .. n], b -> 1))]);
  elif m > 1 and n = 1 then
    return Set([0 .. m], a -> [Concatenation(List([1 .. a], b -> [0]),
                                              List([a + 1 .. m], b -> [1]))]);

  fi;

  row_space := Cartesian(List([1 .. n], i -> [0, 1]));
  row_space := List([0 .. n], i -> Filtered(row_space, row -> Sum(row) = i)); 
  size_row_space := List(row_space, Size);

  bt := function(c)
    local s;
    if Length(c) = rows then
      output(c);
      return;
    fi;
    s := first(c);
    while s <> fail do
      bt(s);
      s := next(s);
    od;
  end;
  
  first := function(c)
    local out;
    out := ShallowCopy(c);
    out[Length(c) + 1] := out[Length(c)];
    return out;
  end;

  next := function(s)
    local out;
    if s[Length(s)] >= size_row_space then
      return fail;
    fi;
    out := ShallowCopy(s);
    out[Length(s)] := out[Length(s)] + 1;
    return out;
  end;

  output := function(c)
    Add(shapes, c);
  end;

  G := SymmetricGroup(cols);
  G := Image(ActionHomomorphism(G, row_space, Permuted));;
  act := function(x, g)
    local xg;
    xg := ShallowCopy(x);
    Apply(xg, a -> a ^ g);
    Sort(xg);
    return xg;
  end;

  shapes := [];
  for i in [1 .. size_row_space] do
    bt([i]);
  od;
  orbs := OrbitsDomain(G, shapes, act);

  out := List(orbs, orb -> List(orb[1], j -> row_space[j]));

  return out;
end;
