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

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (L) = a list (of automorphisms) 
-----
----- OUPUT: A list of signs associated to each automorphism in the 
----- input
-----
----- DESCRIPTION: Given a list of automorphisms, returns a list of 
----- equal size where each automorphism is converted to its sign.
--------------------------------------------------------------------
-------------------------------------------------------------------- 
signOfAuts = method();
signOfAuts(List) := (L) -> (
    apply(#L,i->(
    (-1)^(inversionNumber(L#i))
    ))
    )
 


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (L) = a list (of signs) 
-----
----- OUPUT: True/False
-----
----- DESCRIPTION: Given a list of signs of automorphisms, determines
----- if an odd automorphism is present. 
--------------------------------------------------------------------
-------------------------------------------------------------------- 
hasOddAut = method();
hasOddAut(List) := (L) -> (
    any(L,n -> n<0) 
    )



--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (n,r) = integers 
-----
----- OUPUT: A list of the non-vanishing (n,r)-matroids
-----
----- DESCRIPTION: Given parameters (n,r), returns the list of 
----- those (n,r)-matroids which do not admit odd automorphisms.
--------------------------------------------------------------------
-------------------------------------------------------------------- 
basisCnr = method();
basisCnr(Thing,Thing) := (n,r) -> (
    B := {};
    V = allMatroids(n,r);
    for i from 0 to #V-1 do (if hasOddAut(signOfAuts(aut(V#i))) == false then B = append(B,V#i));
    B
    )









