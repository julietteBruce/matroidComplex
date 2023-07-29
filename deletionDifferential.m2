needsPackage "Matroids"

-- Deletion differential of the matroid complex

noncoloops = method();
noncoloops (Matroid) := (M) -> (
    apply(toList(groundSet M - coloops M), x -> set {x})
)

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
    
    -- Initialize zero column vector {0,...,0} as mutable list
    column := new MutableList from apply(#targetBasis, i -> 0);

    -- Rewrite alternating sum of deletion matroid classes in C_{n-1} 
    ---- to the corresponding sum in QQ^{dim C_{n-1}} wrt the std basis
    sumTerms := deletions(M);
    scan(#sumTerms, i -> (
        scan(#targetBasis, j -> (
            if areIsomorphic(sumTerms#i,targetBasis#i) then column#j = column#j + (-1)^i
        ))
    ));
    toList(column)
)

diffMatrix = method();
diffMatrix (List,List) = (sourceBasis, targetBasis) - > (
    apply((
    ))
)
