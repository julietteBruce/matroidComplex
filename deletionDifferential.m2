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

diffMatrixColumn = method();
diffMatrixColumn (Matroid, List) := (M,targetBasis)-> (
    
    -- Initialize zero column = {0,...,0} as multiple list
    column := new MutableList from apply(#targetBasis, i -> 0);

    -- Compute nonzero terms in columns according to target basis
    terms := deletions(M);
    scan(#terms, i -> (
        scan(#targetBasis, j -> (
            if areIsomorphic(terms#i,targetBasis#i) then  column#j = column#j + (-1)^i
        ))
    ))
)

diffMatrix = method();
diffMatrix (List,List) = (B1,B2) - > (
    apply(L,x->(
        x#
    ))
)


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (L) = a list
-----
----- OUPUT: true or false
-----
----- DESCRIPTION: Given a list L={N,K} of two numbers N and K 
----- representing size of the base set and rank of the matroid
----- this returns true if there are no regular, simple, connected 
----- rank K matroids on [N]. 
-----
--------------------------------------------------------------------
--------------------------------------------------------------------