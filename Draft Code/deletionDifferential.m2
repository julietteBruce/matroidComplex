needsPackage "Matroids"
importFrom("SpechtModule", {"permutationSign"})

--------------------------- withoutOddAut --------------------------
--------------------------------------------------------------------
----- INPUT: Matroid
-----
----- OUTPUT: Boolean
-----
----- DESCRIPTION: Returns true if a matroid admits an odd 
----- automorphism, and false otherwise.
--------------------------------------------------------------------
-------------------------------------------------------------------- 
withoutOddAut = method();
withoutOddAut(Matroid) := (M) -> (
    not any(getIsos(M,M), perm -> permutationSign(perm) == -1)
)

--------------------------- rankedBases --------------------------
--------------------------------------------------------------------
----- INPUT: (Number, Number) = (n,r)
-----
----- OUTPUT: List
-----
----- DESCRIPTION: Given a pair (n,r), r <= n <= 9, return a list of
----- matroids in allMatroids(n,r) without odd automorphisms. This
----- is a basis for C_n^r.
--------------------------------------------------------------------
-------------------------------------------------------------------- 
rankedBases = method();
rankedBases(Number,Number) := (n,r) -> (
    select(allMatroids(n,r), withoutOddAut)
)

-------------------------- diffMatrixColumn ------------------------
--------------------------------------------------------------------
----- INPUT: Matroid
-----
----- OUTPUT: List
-----
----- DESCRIPTION: Given a matroid representing a basis element of (a
----- subspace of) C_n, return the corresponding column vector of the
----- matrix representing the deletion differential from C_n to C_{n-1}
----- with respect to the standard bases on Q^{dim C_n} and 
----- Q^{dim C_{n-1}}.
--------------------------------------------------------------------
-------------------------------------------------------------------- 

diffMatrixColumn = method();
diffMatrixColumn (Matroid, List) := (M,targetBasis)-> (
    
    ---- Initialize zero column vector {0,...,0} as mutable list
    column := new MutableList from apply(#targetBasis, i -> 0);

    noncoloops := toList(groundSet M - coloops M);

    -- Rewrite alternating sum of deletion matroid classes in C_{n-1} 
    -- to the corresponding sum in QQ^{dim C_{n-1}} wrt the std basis
    scan(noncoloops, e -> (
        scan(#targetBasis, j -> (
            if areIsomorphic(deletion(M, set {e}),targetBasis#j) then (
                column#j = column#j + (-1)^(e+1)
                )
        ))
    ));
    -- return immutable list
    toList(column)
)

----------------------------- diffMatrix ---------------------------
--------------------------------------------------------------------
----- INPUT: (List,List)
-----
----- OUTPUT: Matrix
-----
----- DESCRIPTION: Given the pair of a basis for a subspace V of C_n  
----- and a basis for a subspace W of C_{n-1} containing the image
----- of V under the deletion differential, return the matrix 
----- representing the standard bases on QQ^{dim V} and QQ^{dim W},
----- respectively. The first element of the first list will
----- correspond to the basis vector e_1 = (1,0,...,0), etc.
--------------------------------------------------------------------
-------------------------------------------------------------------- 

diffMatrix = method();
diffMatrix (List,List) := (sourceBasis, targetBasis) -> (
    transpose matrix apply(sourceBasis, M -> diffMatrixColumn(M,targetBasis))
)

----------------------------- diffMaps -----------------------------
--------------------------------------------------------------------
----- INPUT: List
-----
----- OUTPUT: List
-----
----- DESCRIPTION: Given an ordered list of matroid bases for a 
----- subcomplex of C_* of the form 
-----    {basis of C_n,  basis of C_{n-1}, ... , basis of C_{n-k}},
----- return the corresponding list of deletion differential matrices
-----    {D_n, D_{n-1}, ... , D_{n-k+1}}.
----- This can then become the input for M2's built in chainComplex
----- method.
--------------------------------------------------------------------
-------------------------------------------------------------------- 

diffMaps = method();
diffMaps (List) := (subcomplexBases)  -> (
    apply( #subcomplexBases - 1, 
        i -> diffMatrix(subcomplexBases#i,subcomplexBases#(i + 1))
    )
)
