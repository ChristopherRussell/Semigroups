# I assume points in 3D space are represented as [x,y,z], where x runs from [1..Dimensions[1]], y runs from [1..Dimensions[2]], and z runs from [1..Dimensions[3]].

# These two functions just map between a 3d matrix and a 1d list, which is a little annoying as GAP counts from 1
Flatten3DPoint := function(Dimensions, point)
   return (point[1]-1) * Dimensions[2] * Dimensions[3] + (point[2]-1) * Dimensions[3] + (point[3]-1) + 1;
end;

Unflatten3DPoint := function(Dimensions, value)
   local ret;
   ret := [];
   value := value - 1;
   ret[3] := value mod Dimensions[3] + 1;
   value := value - (ret[3] - 1);
   value := value / Dimensions[3];
   ret[2] := value mod Dimensions[2] + 1;
   value := value - (ret[2] - 1);
   ret[1] := value / Dimensions[2] + 1;
   return ret;
end;

# Helpers for working between binary matrix shapes and 3d representations
ShapeTo3D := function(shape, Dimensions)
  local 3dshape, point, i;
  3dshape := [];
  for i in [1 .. Dimensions[1] * Dimensions[2]] do
    point := Unflatten2DPointIn3D(Dimensions, i);
    if i in shape then
      point[3] := 2;
      Add(3dshape, point);
    else 
      Add(3dshape, point);
    fi;
  od;
  Apply(3dshape, a -> Flatten3DPoint(dim, a));
  return 3dshape;
end;
  
Flatten3DPointIn2D := function(Dimensions, point)
  return (point[1]-1) * Dimensions[2] + point[2];
end;

Unflatten2DPointIn3D := function(Dimensions, value)
   local ret;
   ret := [];
   value := value - 1;
   ret[2] := value mod Dimensions[2] + 1;
   value := value - (ret[2] - 1);
   value := value / Dimensions[2];
   ret[1] := value mod Dimensions[1] + 1;
   ret[3] := 1;
   return ret;
end;

# Given a 3D cube, dimensions Dimensions, apply permutation 'perm' to dimension 'dim'
ApplyPermWholeDimension := function(Dimensions, dim, perm)
   local map, point, i;
   map := [];
   for i in [1..Product(Dimensions)] do
       point := Unflatten3DPoint(Dimensions, i);
       point[dim] := point[dim]^perm;
       map[i] := Flatten3DPoint(Dimensions, point);
   od;
   return PermList(map);
end;

# Given a 3D cube, dimensions Dimensions, apply permutation 'perm' to dimension 'dim'
# But only to elements which have value 'fixval' in dimension 'fixdim'
ApplyPermSingleAssignDimension := function(Dimensions, dim, perm, fixdim, fixval)
   local map, point, i;
   map := [];
   for i in [1..Product(Dimensions)] do
       point := Unflatten3DPoint(Dimensions, i);
       if point[fixdim] = fixval then
           point[dim] := point[dim]^perm;
       fi;
       map[i] := Flatten3DPoint(Dimensions, point);
   od;
   return PermList(map);
end;

# Assumed that rows and columns are dimensions 1 and 2.
TranspositionPerm := function(Dimensions)
  local map, point, i;
  map := [];
  for i in [1..Product(Dimensions)] do
    point := Unflatten3DPoint(Dimensions, i);
    map[i] := Flatten3DPoint(Dimensions, [point[2], point[1], point[3]]);
  od;
  return PermList(map);
end;

RZMSMatrixIsomorphismGroup := function(shape, nr_rows, nr_cols, G)
  local rows, cols, S, H, gens, elms, rmlt, lmlt, grswaps, gcswaps, g;
  # Row Swaps
  S := SymmetricGroup(nr_rows);
  rows := List(GeneratorsOfGroup(S), x -> ApplyPermWholeDimension(Dim, 1, x));
  
  # Col swaps
  S := SymmetricGroup(nr_cols);  
  cols := List(GeneratorsOfGroup(S), x -> ApplyPermWholeDimension(Dim, 2, x));

  # Subgroup of swaps which stablizes shape
  dim := [nr_rows, nr_cols, Size(G) + 1];
  3dshape := ShapeTo3D(shape, dim);
  H := Stabilizer(Group(Flat([cols, rows])), 3dshape, OnSets);


  gens := GeneratorsOfGroup(G);
  elms := ShallowCopy(Elements(G));

  # Apply g to each row (right mult):
  rmlt := List(gens, g -> PermList(Concatenation([1],
          1 + List(elms, e -> Position(elms, e * g)))));
  grswaps := List([1..Dim[1]], r -> List(rmlt, g ->
  ApplyPermSingleAssignDimension(Dim, 3, g, 1, r)));
  
  # Apply g to each col (left mult by inverse):
  lmlt := List(gens, g -> PermList(Concatenation([1],
          1 + List(elms, e -> Position(elms, g ^ -1 * e)))));
  gcswaps := List([1..Dim[2]], r -> List(lmlt, g ->
  ApplyPermSingleAssignDimension(Dim, 3, g, 2, r)));

  # The RZMS matrix isomorphism group
  g := Group(Flat([GeneratorsOfGroup(H), grswaps, gcswaps]));
  
  #if not Size(g) =  Factorial(dim[1]) * Factorial(dim[2]) * (dim[3] - 1) ^
  #  dim[1] * (dim[3] - 1) ^ dim[2] then
  # Error();
  #fi;

  return g;
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
RZMSMatrices := function(nr_rows, nr_cols, G)
  local out, i, shape, shapes;

  # The m x n case is deduced from the n x m case.
  if nr_rows < nr_cols then
    out := RZMSMatrices(nr_rows, nr_cols, G);
    Apply(out, a -> Unflatten3DPoint([nr_rows, nr_cols, Size(G) + 1], a));
    for i in out do
      i := [i[2], i[1], i[3]];
    od;
    Apply(out, a -> Flatten3DPoint([nr_rows, nr_cols, Size(G) + 1], a));
    return out;
  fi;

  # Get shapes
  shapes := BinaryMatrixShapesByDim(nr_rows, nr_cols);

  # Find by shape
  out := [];
  for shape in shapes do
    Add(out, RZMSMatricesByShape(nr_rows, nr_cols, G, shape));
  od;

  return out;
end;

RZMSMatricesByShape := function(nr_rows, nr_cols, G, shape)
  local dim, space, iter;

  IG := RZMSMatrixIsomorphismGroup(shape, nr_rows, nr_cols, G);
  dim := [nr_rows, nr_cols, Size(G) + 1];
  space := [];
  for i in [1 .. dim[1] * dim[2]] do
    point := Unflatten2DPointIn3D(dim, i);
    if i in shape then
      Add(space, List([1 .. Size(G)], g ->
        Flatten3DPoint(dim, [point[1], point[2], g + 1])));
    else
      Add(space, [Flatten3DPoint(dim, [point[1], point[2], 1])]);
    fi;
  od;

  out := Set([]);
  iter := IteratorOfCartesianProduct(space);
  while not IsDoneIterator(iter) do
    mat := NextIterator(iter);
    AddSet(out, CanonicalImage(IG, mat, OnSets));
  od;

  return out;
end;

