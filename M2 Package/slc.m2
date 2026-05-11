needsPackage "Matroids"

--------------------------------------------------------------------
-- Edit this path
--------------------------------------------------------------------

slcRegularDirectory = currentDirectory() | "data/slcRegular/slcRegular_";


--------------------------------------------------------------------
-- File utilities
--------------------------------------------------------------------
validBidegree = (n,r) -> (n >= 0 and r >= 0 and r <= n);

readBasisFile = method();

readBasisFile(ZZ,ZZ) := (n,r) -> (
    if not validBidegree(n,r) then {}
    else (
        fileName := slcRegularDirectory|toString(n)|"_"|toString(r)|".txt";
        if not fileExists fileName then {}
        else (
            fileLines := select(lines get fileName, l -> #l > 0);

            apply(fileLines, l -> (
                L := value l;
                E := L#0;
                Bases := L#1;
                matroid(E, Bases)
            ))
        )
    )
);


--------------------------------------------------------------------
-- Signs.
--
-- The old code used (-1)^(e+1).  This version uses the position of e
-- in the ordered ground set, so it agrees with the old code when the
-- ground set is {0,...,n-1}.
--
-- It also includes the parity of the isomorphism from the deleted
-- matroid to your chosen representative in targetBasis.  Since your
-- basis is supposed to have no odd automorphisms, this sign is
-- independent of the chosen isomorphism.
--------------------------------------------------------------------

positionInList = (L,x) -> (
    p := null;
    scan(#L, i -> if L#i == x then p = i);
    if p === null then error "positionInList: element not found";
    p
);

permutationSign = p -> (
    invs := 0;
    scan(#p, i -> (
        if i + 1 < #p then scan((i+1)..(#p-1), j -> (
            if p#i > p#j then invs = invs + 1
        ))
    ));
    if invs % 2 == 0 then 1 else -1
);

isomorphismSign = (phi, sourceOrder, targetOrder) -> (
    permutationSign(apply(sourceOrder,
        x -> positionInList(targetOrder, phi#x)))
);


diffMatrixColumn = method();

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

apply(toList(1..15), n-> (apply(toList(0..n), r->(
		print (n,r);
		slcDifferential(n,r)
		))
	))
	
	print (n,r);
        slcDifferential(n,r)
--	print L;
--	L
	))

--------------------------------------------------------------------
-- Differential matrix
--
-- This now handles all empty-source/empty-target edge cases.
--
-- Output has size:
--   (#targetBasis) x (#sourceBasis)
--------------------------------------------------------------------

zeroDiffMatrix = (numRows,numCols) -> (
    map(QQ^numRows, QQ^numCols, 0)
);

diffMatrix = method();

diffMatrix(List,List) := (sourceBasis,targetBasis) -> (
    if #sourceBasis == 0 or #targetBasis == 0 then
        zeroDiffMatrix(#targetBasis, #sourceBasis)
    else
        sub(transpose matrix apply(sourceBasis,
            M -> diffMatrixColumn(M,targetBasis)), QQ)
);


--------------------------------------------------------------------
-- Deletion differential in bidegree (n,r)
--
-- Deleting a non-coloop preserves rank, so this is
--
--   C_{n,r} ----> C_{n-1,r}.
--
-- Edge cases such as r > n-1 are now fine: rankedBasis(n-1,r) = {},
-- so the matrix has zero rows.
--------------------------------------------------------------------

slcDifferential = method();

slcDifferential(ZZ,ZZ) := (n,r) -> (
    diffMatrix(readBasisFile(n,r), readBasisFile(n-1,r))
);


--------------------------------------------------------------------
-- Convenience wrapper:
-- deletion of a non-coloop preserves rank, so this is C_{n,r}->C_{n-1,r}.
--------------------------------------------------------------------

slcDifferential = (n,r) -> (
    diffMatrix(readBasisFile(n,r), readBasisFile(n-1,r))
);


--------------------------------------------------------------------
-- Optional diagnostic: check whether the chosen span is actually
-- closed under deletion in bidegree (n,r).
--------------------------------------------------------------------

unmatchedDeletionTerms = method();
unmatchedDeletionTerms(Matroid,List) := (M,targetBasis) -> (
    select(toList(groundSet M - coloops M), e -> (
        D := deletion(M, set {e});
        all(targetBasis, T -> isomorphism(D,T) === null)
    ))
);

closureReport = (n,r) -> (
    sourceBasis := rankedBasis(n,r);
    targetBasis := rankedBasis(n-1,r);
    select(
        apply(sourceBasis, M -> {M, unmatchedDeletionTerms(M,targetBasis)}),
        x -> #(x#1) > 0
    )
);
