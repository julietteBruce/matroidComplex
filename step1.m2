needsPackage "Matroids";

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (M) = a matroid 
-----
----- OUPUT: A list of the automorphisms of the matroid.
-----
----- DESCRIPTION: Given a matroid , a list give the Fripertinger-Wild
----- representation of a representing a matroid, or a matrix
----- representing a matroid this function returns a list of the 
----- automorphisms of the matroid.
----- 
--------------------------------------------------------------------
--------------------------------------------------------------------

aut = method();
aut(Thing) := (M) -> (
    getIsos(M,M)
    )

---- Test ----
F7 = specificMatroid "fano"
assert(#aut(F7) == 168)
--------------

isIsomorphicMatroids = method();
isIsomorphicMatroids(Thing, Thing) := (M,N) -> (
    A = getIsos(M,N);
    A != {}
    )

inversionNumber = method();
inversionNumber(List) := (P) -> (
    sum flatten apply(#P,i->(
	 delete(,apply(toList(i..#P-1),j->(
		if P#j < P#i then 1
		)))
	))
    )

P = {3,2,1}
inversionNumber(P)
PP = {0,6,5,4,3,2,1}
inversionNumber(PP)

sign = method();
sign(List) := (P) -> (
    (-1)^(inversionNumber(P))
    )

P = {3,2,1}

sum flatten apply(#P,i->(
	 delete(,apply(toList(i..#P-1),j->(
		if P#j < P#i then 1
		)))
	))
