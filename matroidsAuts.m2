needsPackage "Matroids";
needsPackage "TestIdeals";
needsPackage "SpechtModule";
needsPackage "TensorComplexes";



--------------------------------------------------------------------
--------------------------------------------------------------------
----- NOTES
-----
----- #0 The Fripertinger-Wild database seemed to be missing the 
----- file for N=10, K=9. I guessed what this should be and added it.
-----
----- #1 The Fripertinger-Wild database only has simple, connected,
----- regular, LOOPLESS, matroids and we might need non-loopless
----- matroids??? 
----- Dan: No, the "non-loopless" is redundant, simple matroids have 
----- no loops.
--------------------------------------------------------------------
--------------------------------------------------------------------



--------------------------------------------------------------------
--------------------------------------------------------------------
----- DESCRIPTION: This snippet of code goes through all of the files
----- in the simpleRegularLooplessMatroids folder, and extracts the 
----- size of the base set and the rank for each given file. For
----- example, reg_simple_con3_2.txt is converted to the list {3,2}.
----- The output is a list of these pairs called dataRange.
--------------------------------------------------------------------
--------------------------------------------------------------------
fileList = lines get "! ls simpleRegularLooplessMatroids";
dataRange = apply(fileList,s->(
	S := separateRegexp("[_]",s);
	L := separateRegexp("[.]",S#3);
	{value substring(3,(S#2)),value L#0}
	));



--------------------------------------------------------------------
--------------------------------------------------------------------
----- DESCRIPTION: This is a list of pairs {N,K} for K<=N-1 where we 
----- know that there are no simple, connected, regular matroids on 
----- [N] of rank K. There are no such matroids when K>=N.
--------------------------------------------------------------------
--------------------------------------------------------------------
dataZero = {{2,1},{3,1},{4,1},{4,2},{5,1},{5,2},{6,1},{6,2},{7,1},{7,2},
    {7,3},{8,1},{8,2},{8,3},{9,1},{9,2},{9,3},{10,1},{10,2},{10,3},{11,1},
    {11,2},{11,3},{11,4},{12,1},{12,2},{12,3},{12,4},{13,1},{13,2},{13,3},
    {13,4},{14,1},{14,2},{14,3},{14,4},{15,1},{15,2},{15,3},{15,4}};



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

noMatroidsNK = method();
noMatroidsNK (List) := (L) ->(
    N:=L#0; K:=L#1; 
    if (N>=2 and K>=N) then return true;
    if (N==1 and K>N) then return true;
    if member(L,dataZero) then return true;
    if not member(L,dataRange) then return "no data for this {N,K}";
    return false
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (L) = a list
-----
----- OUPUT: A string representing a file name
-----
----- DESCRIPTION: Given a list L={N,K} of two numbers N and K 
----- representing size of the base set  and rank of the matroid
----- this returns the file name and path for the corresponding 
----- file in the simpleRegularLooplessMatroids folder as a string. 
----- Specifically it returns
----- "simpleRegularLooplessMatroids/reg_simple_conN_K" 
----- where N and K have been converted to strings. 
-----
--------------------------------------------------------------------
--------------------------------------------------------------------

dataFileNameFromList = method();
dataFileNameFromList (List) := (L) ->(
    "simpleRegularLooplessMatroids/reg_simple_con"|toString(L#0)|"_"|toString(L#1)|".txt"
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (L) = a list representing a matroid
-----
----- OUPUT: A matrix representing the matroid
-----
----- DESCRIPTION: Given a list L = {L1,L2,...,Ln} of positive intgers
----- which represents a matroid in Fripertinger-Wild form, this
----- returns the matrix M whose i-th column is the base 2 expansion
----- of Li written increasingly and padded with zeros so every 
----- element of L has a base 2 expansion of the same length. The
----- resulting matrix M represents the matroid over GF(2).
-----
--------------------------------------------------------------------
--------------------------------------------------------------------

convertListToMatrix = method();
convertListToMatrix (List) := (L) ->(
    L1 := apply(L,i->(
	    adicExpansion(2,i)
	    ));
    maxLen := max apply(L1,i->#i);
    L2 := apply(L1,i->(
	    i|apply((maxLen-#i),j->0)
	    ));
    return sub(transpose matrix L2, GF(2))
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (L) = a list representing a matroid
-----
----- OUPUT: A matroid
-----
----- DESCRIPTION: Given a list L = {L1,L2,...,Ln} of positive intgers
----- which represents a matroid in Fripertinger-Wild form, this
----- returns matroid corresponding to the matrix from the function 
----- "convertListToMatrix" 
-----
--------------------------------------------------------------------
--------------------------------------------------------------------

convertListToMatroid = method();
convertListToMatroid (List) := (L) -> matroid(convertListToMatrix(L))

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (M) = a matrix representing a matroid
-----
----- OUPUT: A list representing the matroid
-----
----- DESCRIPTION: Given a binary matrix M with columns (M1,M2,...Mn) 
----- representing a matroid over RR this returns a list of positive
----- integers L = {L1,L2,...,Ln} such that the i-th  column is the
----- base 2 expansion of Li written increasingly and padded with 
----- zeros so every element of L has a base 2 expansion of the same 
----- length. The list L is the Fripertinger-Wild representation of
----- of the matroid represented by M.
-----
--------------------------------------------------------------------
--------------------------------------------------------------------

convertMatrixToList = method();
convertMatrixToList (Matrix) := (M) ->(
    L1 := entries transpose M;
    apply(L1,N->(
	    sum(#N,j->(N#j)*2^(j))
	    ))
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (L) = a list 
-----
----- OUPUT: A list of matroids
-----
----- DESCRIPTION: Given a list L={N,K} this function returns a list
----- of all the regular, simple, connected matroids 
----- whose base set has size N with rank K,  up to isomorphism. 
----- This is done by loading the matroids from the corresponding 
----- data file generated from the Fripertinger-Wild data - in the 
----- simpleRegularLooplessMatroids folder. 
-----
--------------------------------------------------------------------
--------------------------------------------------------------------

rslcMatroids = method();
rslcMatroids (List) := (L) ->(
    if noMatroidsNK(L) then return {};
    if member(L,dataRange) then (
	L1 := value get dataFileNameFromList(L);
    	return apply(L1,i->convertListToMatroid(i) );
	)
    else return "No data for this {N,K}."
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: () = a positive integer
-----
----- OUPUT: A list of matroids
-----
----- DESCRIPTION: Given a positive integer N, returns all regular, 
----- simple connected matroids on [N]. 
-----
--------------------------------------------------------------------
--------------------------------------------------------------------

rslcMatroidsOnN = method();
rslcMatroidsOnN (ZZ) := (N) -> (
    if N<=0 then return "N is not positive";
    if N>=16 then return "no data";
    return flatten for K in (1..N) list rslcMatroids({N,K});
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: ZZ,ZZ
-----
----- OUPUT: a List 
-----
----- DESCRIPTION: given integers a,b, generates all size a multisets  
----- of regular, simple, connected matroids on [b].
--------------------------------------------------------------------
--------------------------------------------------------------------

partitionsSameGroundSet = method();
partitionsSameGroundSet (ZZ,ZZ) := (a,b) -> (    
    allMatsOnb := rslcMatroidsOnN(b);
    return multiSubsets(allMatsOnb, a)
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: List
-----
----- OUPUT: List 
-----
----- DESCRIPTION: given a list of sets, computes the cartesian 
----- product of them 
----- 
--------------------------------------------------------------------
--------------------------------------------------------------------

cartesianProductSets = method();
cartesianProductSets (List) := (L) -> (    
    if #L==1 then return apply(toList(L#0), j->{j});    
    A := L#0;
    for i in (1..#L-1) do (
	A=(A ** L#i)/splice;
	);
    return apply(toList(A), j->toList(j))
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: Partition
-----
----- OUPUT: List 
-----
----- DESCRIPTION: given a Partition P={P1,...,Pa} (P1<=...<=Pa), 
----- generates all multisets {M1,...,Ma} of regular, simple,
----- connected matroids Mi on [Pi]. 
----- 
--------------------------------------------------------------------
--------------------------------------------------------------------

matroidsOfPartition = method();
matroidsOfPartition (Partition) := (P) -> (    
    lP := toList(P);
    ulP := unique lP;
    numPart := hashTable(for i in ulP list {i, number(lP, j->i==j)});
    mats := apply(ulP, b->set(partitionsSameGroundSet(numPart#b,b)));
    return apply(cartesianProductSets(mats), J-> flatten J) 
    ) 



--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (L) = a list of matroids
-----
----- OUPUT: A matroid
-----
----- DESCRIPTION: Given a list of matroids, returns their direct sum.
-----
--------------------------------------------------------------------
--------------------------------------------------------------------

directSumMatroids = method();
directSumMatroids (List) := (L) -> (
    if #L==1 then return L#0;
    A := L#0;
    for i in (1..#L-1) do (
	A = A ++ L#i;
	);
    return A
    )    


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (L) = an integer
-----
----- OUPUT: A list of matroids
-----
----- DESCRIPTION: Given an integer this function returns a list
----- of all the regular, simple, matroids on [N], up to 
----- isomorphism. This is done by taking direct sums of smaller 
----- regular, simple, connected matroids.
-----
--------------------------------------------------------------------
--------------------------------------------------------------------

rslMatroidsOnN = method();
rslMatroidsOnN (ZZ) := (N) -> (
    pN := partitions N;
    listMats := flatten apply(pN, P-> matroidsOfPartition(P)); 
    return apply(listMats, Ms -> directSumMatroids(Ms))
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (L) = a list 
-----
----- OUPUT: A list of matroids
-----
----- DESCRIPTION: Given a list L={N,K} this function returns a list
----- of all the regular, simple, matroids whose base set
----- has size N with rank K, up to isomorphism. This is done by 
----- taking direct sums of smaller regular, simple matroids 
-----
--------------------------------------------------------------------
--------------------------------------------------------------------

rslMatroids = method();
rslMatroids (List) := (L) -> (
    rslMonN := rslMatroidsOnN(L#0);
    return for i in (0..#rslMonN-1) list if rank(rslMonN#i)==L#1 then rslMonN#i else continue;
    )


    
--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (M) = a matroid 
----- INPUT: (L) = a list representing a matroid
----- INPUT: (M) = a matrix 
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
aut (Thing) := (M) ->(
    getIsos(M,M)
    )

aut (List) := (L) ->(
    aut(matroid(convertMatrixToList(L)))
    )

aut (Matrix) := (M) ->(
    aut(matroid(M))
    )

---- Test ----
F7 = specificMatroid "fano"
assert(#aut(F7) == 168)
--------------

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (L) = a list
-----
----- OUPUT: A list of the size of automorphism of matroids
-----
----- DESCRIPTION: Given a list L={N,K} this function returns a list
----- chose elements correspond to the size of the automorphism group
----- of all the regular, simple, connected, loopless matroids whose 
----- base set has size N with rank K, up to isomorphism. 
-----
--------------------------------------------------------------------
--------------------------------------------------------------------
rslMatroidsNumAut = method();
rslMatroidsNumAut (List) := (L) ->(
    L1 := rslMatroids(L);
    apply(L1,M->(#aut(M)))
    )

---- seems to bog down around {11,9}
H = new MutableHashTable;
apply(sort dataRange,D->(
	H#D = rslMatroidsNumAut(D);
	print D
	))


rslcAlternatingMatroids = method();
rslcAlternatingMatroids (List) := (L) ->(
    L1 := rslcMatroids(L);
    delete(,apply(L1,M->(
		G := aut(M);
		L2 := unique apply(G,C->(
			if permutationSign(C) == -1 then break {false}
			else true
			));
		if L2 == {true} then M
		)))
       )

apply(sort dataRange, D->(
	L1 := rslcAlternatingMatroids(D);
	print {D,#L1}
	))


{{1, 1}, 1}
{{3, 2}, 0}
{{4, 3}, 0}
{{5, 3}, 0}
{{5, 4}, 0}
{{6, 3}, 1}
{{6, 4}, 0}
{{6, 5}, 0}
{{7, 4}, 0}
{{7, 5}, 0}
{{7, 6}, 0}
{{8, 4}, 0}
{{8, 5}, 0}
{{8, 6}, 0}
{{8, 7}, 0}
{{9, 4}, 0}
{{9, 5}, 0}
{{9, 6}, 0}
{{9, 7}, 0}
{{9, 8}, 0}
{{10, 4}, 0}
{{10, 5}, 2}
{{10, 6}, 0}
{{10, 7}, 0}
{{10, 8}, 0}
{{10, 9}, 0}
{{11, 5}, 1}
{{11, 6}, 1}
{{11, 7}, 0}
{{11, 8}, 0}
{{11, 9}, 0}
{{12, 5}, 0}
{{12, 6}, 5}
{{12, 7}, 0}
{{13, 5}, 0}
{{13, 6}, 12}
{{13, 6}, 12}

rslAlternatingMatroids = method();
rslAlternatingMatroids (List) := (L) ->(
    L1 := rslMatroids(L);
    delete(,apply(L1,M->(
		G := aut(M);
		L2 := unique apply(G,C->(
			if permutationSign(C) == -1 then break {false}
			else true
			));
		if L2 == {true} then M
		)))
       )


rslCoLooplessMatroids = method();
rslCoLooplessMatroids (List) := (L) ->(
    L1 := rslMatroids(L);
    delete(,apply(L1,M->(
		L2 := coloops(M);
		if L2 == {} then M
		)))
    )



matroidChi =  method();
matroidChi (ZZ) := (g) ->(
    L1 := sort delete(,apply(dataRange,D->(
		if (D#1) == g then D
		)));
    sum apply(L1,D->(
	    L2 := rslCoLooplessMatroids(D);
	    L3 := delete(,apply(L2,M->(
			G := aut(M);
			L4 := unique apply(G,C->(
				if permutationSign(C) == -1 then break {false}
				else true
				));
			if L4 == {true} then M
			)));
	    (-1)^(D#0)*(#L3)
	    ))
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (g) = a number
-----
----- OUPUT: "The" orbitfold Euler characterisitc of R_g
-----
----- Description: Given a positive integer g this function returns
----- "the" orbitfold Euler characteristic of R_g i.e. it sums
----- (-1)^(|B(M)|)*1/|Aut(M)| where B(M) is the base set of M such 
----- where the sum is over all connected, simple, loopless, matroids
----- of rank at most g.
----- 
--------------------------------------------------------------------
--------------------------------------------------------------------
matroidOrbChi =  method();
matroidOrbChi (ZZ) := (g) ->(
    L1 := delete(,apply(dataRange,D->(
		if (D#1) <= g then D
		)));
    sum(L1,D->(
	    (-1)^(D#0)*(sum rslMatroidsNumAut(D))^(-1)
	    ))
    )


matroidOrbChi(3)
matroidOrbChi(4)
matroidOrbChi(5)
matroidOrbChi(6)
matroidOrbChi(7)
matroidOrbChi(8)

--- n







--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (M) = a matroid with no loops or coloops
-----
----- OUPUT: a graph; The basis graph of M
-----
----- Description: 
----- 
--------------------------------------------------------------------
--------------------------------------------------------------------
basisGraph =  method();
basisGraph (Thing) := (M) ->(
    if #(loops M)>0 or #(coloops M) >0 then return "Matroid has loops or coloops.";
    V := bases(M); --Vd := hashTable for i in (0..#V-1) list {V#i,i};
    V2 := subsets(#V,2); E := {}; r=rank(M);
    for k in (0..#V2-1) do (
        e := V2#k; e0:=e#0; e1:=e#1;
	if #(V#e0 * V#e1) == r-1 then E=append(E,{V#e0, V#e1 });
	); 
    return graph(E)
    )
	











flatsRk = (M,k) -> (
    F:=flats(M);
    return for i in (0..#F-1) list if rank(M,F#i)==k then F#i else continue
    )


-- given matroid M and file name fName, makes a file named fname that contains the bases of M. 
-- Each line contains a basis as a space separated list.
matroidToBasesFile = (M,fName) -> (
    Bs := apply(bases(M), b->sort toList b);
    f := openOut(fName);
    for b in Bs do (
	bs := apply(b, i->toString(i));
	f << demark(" ", bs) << endl;	
	);
    close f;
    return fName
    )

-- given matroid M and file name fName, makes a file named fname that contains the rank 2 flats of M. 
-- Each line contains a flat as a space separated list.
matroidToRk2FlatsFile = (M,fName) -> (
    Fs := apply(flatsRk(M,2), b->sort toList b);
    f := openOut(fName);
    for g in Fs do (
	gs := apply(g, i->toString(i));
	f << demark(" ", gs) << endl;	
	);
    close f;
    return fName
    )


    
M615s = rslcMatroids({15,6})
for i in (0..#M615s-1) do (matroidToBasesFile(M614s#i, "srlmRank6n15Bases/M"|i|".dat"))
for i in (0..#M615s-1) do (matroidToRk2FlatsFile(M614s#i, "srlmRank6n15Rank2Flats/M"|i|".dat"))

    

-- #(nonbases(M514))
-- time subsets(864,2);
-- binomial(864,2)
-- help Matroid.rank


A=rslcAlternatingMatroids({13,6})
apply(A, M-> #(aut(M)))
