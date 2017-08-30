###############################################################################
# Permutations and operations on sets represenating (regular) binary matrices
###############################################################################

# Calculates the perm group which acts on these sets (via OnSets) by
# transposition, row permutations and column permutations
MatRowPermGrp := function(nr_rows, nr_cols)
  local out, temp, a, i;

  if nr_rows = 1 and nr_cols = 1 then
    return Group(PermList([1]));
  fi;

  out := [];
  if nr_rows = nr_cols then
  temp := [];
  for i in [1 .. nr_rows] do
    Append(temp, List([1 .. nr_cols], j -> i + (j - 1) * nr_cols));
  od;
  Add(out, PermList(temp));
  fi;

  if nr_rows > 1 then
    Add(out, PermList(
      Concatenation([nr_cols + 1 .. 2 * nr_cols], [1 .. nr_cols])));
    if nr_rows > 2 then
      Add(out, PermList(
        Concatenation([(nr_rows - 1) * nr_cols + 1 .. nr_rows * nr_cols],
                      [1 .. (nr_rows - 1) * nr_cols])));
    fi;
  fi;

  if nr_cols > 1 then
    temp := Concatenation([2, 1], [3 .. nr_cols]); 
    Add(out, PermList(
      Concatenation(temp + nr_cols, temp, 
        Concatenation(List([3 .. nr_rows], a -> temp + nr_cols * (a - 1)))))); 
    if nr_cols > 2 then
      temp := Concatenation([2 .. nr_cols], [1]);
      Add(out, PermList(
        Concatenation(List([1 .. nr_rows], a -> temp + nr_cols * (a - 1)))));  
    fi;
  fi;

  return Group(out);
end;

# This permutation sends (via OnSets) a set representing an m x n matrix to its
# transpose as a set representing a n x m matrix.
TranspositionPermutation := function(nr_rows, nr_cols)
  return PermList(Concatenation(
           List([1 .. nr_rows], i ->
             List([1 .. nr_cols], j -> (j - 1) * nr_rows) + i)));
end;

# This returns the binary matrix represented by the set <set>. Note this
# depends on the dimensions!
SetToBinaryMat := function(set, nr_rows, nr_cols, G)
  local func, out, i, row;
  func := function(b)
    if b then
      return 1;
    fi;
    return 0;
  end;
  out := EmptyPlist(nr_rows);
  for i in [1 .. nr_rows] do
    row := List([(i - 1) * nr_cols + 1 .. i * nr_cols], a -> func(a in set));
    Add(out, ShallowCopy(row));
  od;
  return out;
end;

# This returns a transformation which relabels the elements of the set so that
# it represents the upper left corner of a matrix with more rows and/or columns
EmbedSetTransformation := function(nr_rows, nr_cols, new_rows, new_cols)
  local trans;
  trans := Concatenation(List([1 .. nr_rows], i -> List([1 .. nr_cols], j ->
           j + (i - 1) * new_cols)));
  Append(trans, List([1 .. (new_cols - nr_cols) * (nr_rows - 1)], a -> 1));
  return TransformationList(trans);
end;

# Returns the extensions for a n-1 x n-1 01matrix (as a set) to a n x n matrix.
SquareExtensions := function(n)
  local last_col, last_row, ext1, ext2, comb1, comb2, c;
  last_col := Set([1 .. n - 1], i -> i * n);
  last_row := Set([1 .. n - 1], i -> n * (n - 1) + i);

  ext1 := Combinations(Union(last_col, last_row));
  for c in ext1 do
    AddSet(c, n ^ 2);
  od;

  comb1 := Combinations(last_col);
  Remove(comb1, 1);
  comb2 := Combinations(last_row);
  Remove(comb2, 1);

  ext2 := Cartesian(comb1, comb2);
  Apply(ext2, Union);

  return Union(ext1, ext2);
end;

# Returns the extensions for a m-1 x n 01matirx (as a set) to a m x n matrix.
RowExtensions := function(nr_rows, nr_cols)
  local last_row, comb, c;
  last_row := Set([1 .. nr_cols], i -> nr_cols * (nr_rows - 1) + i);

  comb := Combinations(last_row);
  Remove(comb, 1); # Remove the empty row.

  return comb;
end;


###############################################################################
# Finding the binary matrix shapes
###############################################################################

# Solution to one row & one column case.
BinaryMatrixShapesWithOneRowOrCol := function(dim)
  return Set([1 .. dim], a -> Set([1 .. a]));
end;

# Solution to two row case.
BinaryMatrixShapesWithTwoRows := function(nr_cols)
  local row1, row2, out, i, j, opp;
  out := [];
  row1 := [1 .. nr_cols];
  for j in [1 .. nr_cols] do
    Add(out, Union(row1, [1 .. j] + nr_cols));
  od;

  i := nr_cols - 1;
  while 2 * i >= nr_cols do
    row1 := [1 .. i];
    opp := [1 .. nr_cols];
    SubtractSet(opp, row1);
    opp := opp + nr_cols;
    for j in [0 .. 2 * i - nr_cols] do
      row2 := Union(opp, [1 .. j] + nr_cols);
      Add(out, Union(row1, row2));
    od;
    i := i - 1;
  od;
  return out;
end;

# Wrapper for solved cases.
BinaryMatrixShapesSolvedCases := function(nr_rows, nr_cols)
  local out, transpose;
  if nr_rows = 1 or nr_cols = 1 then
    return BinaryMatrixShapesWithOneRowOrCol(Maximum(nr_rows, nr_cols));
  fi;

  if nr_rows = 2 then
    return BinaryMatrixShapesWithTwoRows(nr_cols);
  fi;
end;

###############################################################################
# IO methods used whilst seraching as well as to save results.
###############################################################################
BinaryMatrixShapesReadFile := function(nr_rows, nr_cols)
  local prefix, name, file, out;
  prefix := "/Users/crussell/Desktop/shapes/shapes-";
  name := Concatenation(prefix, String(nr_rows), "x", String(nr_cols));
  file := IO_File(name, "r");
  if file = fail then
    return fail;
  fi;
  out := IO_ReadLines(file);
  IO_Close(file);
  Apply(out, a -> EvalString(Chomp(a)));  
  return out;
end;

BinaryMatrixShapesWriteFile := function(nr_rows, nr_cols, data)
  local prefix, name, file, d;
  prefix := "/Users/crussell/Desktop/shapes/shapes-";
  name := Concatenation(prefix, String(nr_rows), "x", String(nr_cols));
  file := IO_File(name, "w");
  for d in data do
    IO_WriteLine(file, d);
  od;
  IO_Close(file);
end;

###############################################################################
# Main method:
# This finds shapes by building up from lower dimension cases.
###############################################################################
BinaryMatrixShapesByDim := function(nr_rows, nr_cols)
  local out, trans, embed, data, extensions, G, d, e, temp;

  # If we have calculated this case before then return those results
  data := BinaryMatrixShapesReadFile(nr_rows, nr_cols);
  if data <> fail then
    return data;
  fi;

  # These are cases which we can describe the solution to explicitly
  if nr_rows <= 2 or nr_cols = 1 then
    return BinaryMatrixShapesSolvedCases(nr_rows, nr_cols);
  fi;

  # The m x n case is deduced from the n x m case.
  if nr_rows < nr_cols then
    out := BinaryMatrixShapesByDim(nr_cols, nr_rows);
    trans := TranspositionPermutation(nr_rows, nr_cols);
    Apply(out, a -> OnSets(a, trans));

  # In square matrix case, build up from previous square case
  elif nr_rows = nr_cols then
    data := BinaryMatrixShapesByDim(nr_rows - 1, nr_cols - 1);
    embed := EmbedSetTransformation(nr_rows - 1, nr_cols - 1, nr_rows, nr_cols);
    Apply(data, d -> OnSets(d, embed));
    G := MatRowPermGrp(nr_rows, nr_cols);
    extensions := SquareExtensions(nr_rows);
    out := [];
    for d in data do
      for e in extensions do
        temp := ShallowCopy(d);
        UniteSet(temp, e);
        AddSet(out, CanonicalImage(G, temp, OnSets));
      od;
    od;

  # If more rows than columns then build from a square matrix by adding rows.
  elif nr_rows > nr_cols then
    data := BinaryMatrixShapesByDim(nr_rows - 1, nr_cols);
    G := MatRowPermGrp(nr_rows, nr_cols);
    extensions := RowExtensions(nr_rows, nr_cols);
    out := EmptyPlist(Size(extensions) * Size(data));
    for d in data do
      for e in extensions do
        temp := ShallowCopy(d);
        UniteSet(temp, e);
        AddSet(out, CanonicalImage(G, temp, OnSets));
      od;
    od;
  fi;

  BinaryMatrixShapesWriteFile(nr_rows, nr_cols, out);  
  return out;
end;

###############################################################################
# Viewing matrices:
###############################################################################
PrintMat := function(m)
  for r in m do
    Print(r, "\n");
  od;
  Print("\n");
end;

###############################################################################
# Notes:
###############################################################################
# -Use CanonicalImage with OnSets. [[1,0],[1,1]] = {1,3,4} etc.
# -Build row by row.
# -Use tricks to reduce search space, such as rows need to have less 1's that above
# row. Only do this where worthwhile!
# -
# -Could use CanonicalDigraph with m * n two coloured clique graph which CJ showed

