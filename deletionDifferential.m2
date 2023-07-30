needsPackage "Matroids"

-------------------------- diffMatrixColumn ------------------------
--------------------------------------------------------------------
----- INPUT: Matroid
-----
----- OUTPUT: List
-----
----- DESCRIPTION: Given a Matroid, return a list of coloops as 
----- singletons.
--------------------------------------------------------------------
-------------------------------------------------------------------- 

noncoloops = method();
noncoloops (Matroid) := (M) -> (
    apply(toList(groundSet M - coloops M), x -> set {x})
)

-------------------------- diffMatrixColumn ------------------------
--------------------------------------------------------------------
----- INPUT: Matroid
-----
----- OUTPUT: List
-----
----- DESCRIPTION: Given a Matroid M, return a list of deletion 
----- matroids corresponding to each nonloop in the ground set of M.
--------------------------------------------------------------------
-------------------------------------------------------------------- 

deletions = method();
deletions (Matroid) := (M) -> (
    apply(noncoloops(M), e -> deletion(M,e))
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

    -- Rewrite alternating sum of deletion matroid classes in C_{n-1} 
    -- to the corresponding sum in QQ^{dim C_{n-1}} wrt the std basis
    sumTerms := deletions(M);
    scan(#sumTerms, i -> (
        scan(#targetBasis, j -> (
            if areIsomorphic(sumTerms#i,targetBasis#i) then (
                column#j = column#j + (-1)^i
                )
        ))
    ));
    -- return unmutiple list
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
diffMatrix (List,List) = (sourceBasis, targetBasis) - > (
    transpose matrix apply(sourceBasis, M -> diffMatrixColumn(M,targetBasis))
)
