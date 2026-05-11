-- -*- coding: utf-8 -*-
--------------------------------------------------------------------------------
-- Copyright 2021  Juliette Bruce, Benjamin Ashlock, Jacob Bucciarelli, Bailee Zacovic
--
-- This program is free software: you can redistribute it and/or modify it under
-- the terms of the GNU General Public License as published by the Free Software
-- Foundation, either version 3 of the License, or (at your option) any later
-- version.
--
-- This program is distributed in the hope that it will be useful, but WITHOUT
-- ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
-- FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
-- details.
--
-- You should have received a copy of the GNU General Public License along with
-- this program.  If not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- PURPOSE : Tools for computing the matroid chain complex
--
--
-- PROGRAMMERS : Juliette Bruce, Benjamin Ashlock, Jacob Bucciarelli, Bailee Zacovic
--
--
-- UPDATE HISTORY #0 - 
--
--
-- UPDATE HISTORY #1 - August 2023 - Bailee Zacovic: Began preparing
-- package for eventual publication. Adding tests, comments, documentation,
-- cleaning up code, etc.
--
--
-- UPDATE HISTORY #2 - 
--
--
-- TO DO LIST : create tests
--------------------------------------------------------------------------------



newPackage("MatroidComplexes",
    Version => "1.0",
    Date => "01 August 2023",
    Headline => "Tools for computing the matroid chain complex",
    Authors => {
        {
            Name => "Juliette Bruce",
            Email => "juliette.bruce@berkeley.edu",
            HomePage => "https://juliettebruce.github.io"
        },
        {
            Name => "Benjamin Ashlock",
            Email => "bak6t@missouri.edu"
        },
        {
            Name => "Jacob Bucciarelli",
            Email => "jbucciarelli@ksu.edu"
        },	     
        {
            Name => "Bailee Zacovic",          
            Email => "bzacovic@umich.edu"
    }},
  PackageExports => {"Matroids","SpechtModule"},
  DebuggingMode => true,
  AuxiliaryFiles => true
  )

export {
  "withoutOddAut", 
  "rankedBasis", 
  "diffMatrixColumn", 
  "diffMatrix", 
  "signHashPerm",
  "homologyRank",
  "createBasisFile",
  "readBasisFile",
  "createSubBasisFile",
  "isGraphic",
  "isCographic",
  "isRegularMatroid",
  "isTernary",
  "createDiffMatrixFile"
  }

--------------------------------------------------------------------
--------------------------------------------------------------------
----- CODE
--------------------------------------------------------------------
--------------------------------------------------------------------
-*
needsPackage "Matroids"
importFrom("SpechtModule", {"permutationSign"})
*-

--------------- withoutOddAutsignHashPerm --------------------------
--------------------------------------------------------------------

----- INPUT: Hash Table Representing a Permutation
-----
----- OUTPUT: +1 or -1
-----
----- DESCRIPTION: Returns sign of a permutation stored as a hash.
--------------------------------------------------------------------
-------------------------------------------------------------------- 
signHashPerm = p -> (
    dom := sort keys p;
    imgs := apply(dom, i -> p#i);
    --
    inv := 0;
    n := #imgs;
    --
    for i from 0 to n-2 do (
        for j from i+1 to n-1 do (
            if imgs#i > imgs#j then inv = inv + 1
        )
    );
    --
    if even inv then 1 else -1
)


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

--------------------------- rankedBasis --------------------------
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
rankedBasis = method();
rankedBasis(ZZ, ZZ) := (n, r) -> (
    if r > n then error "Rank must be smaller than ground set.";
    select(allMatroids(n,r), withoutOddAut)
)

-------------------------- diffMatrixColumn ------------------------
--------------------------------------------------------------------
----- INPUT: Matroid, List
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
    --
    noncoloops := toList(groundSet M - coloops M);
    -- Rewrite alternating sum of deletion matroid classes in C_{n-1} 
    -- to the corresponding sum in QQ^{dim C_{n-1}} wrt the std basis
    scan(noncoloops, e -> (
        D := deletion(M, set {e});
        scan(#targetBasis, j -> (
            isos := isomorphism(D, targetBasis#j);
            if isos =!= null then (
                -- Correct paper convention for Macaulay2's 0-based labels:
                column#j = column#j + (-1)^e * signHashPerm(isos)
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

----------------------------- diffMatrix ---------------------------
--------------------------------------------------------------------
----- INPUT: (n,r,string)
-----
----- OUTPUT: Matrix
-----
----- DESCRIPTION: Given integers n, r, and string s this will
----- comput the differential from s_{n,r} ---> s_{n-1,r}
----- where s is the subcomplex referened by the property s.
-----
----- Note s must be one of the strings below:
----- "basis"
----- "binary"
----- "cographic"
----- "graphic"
----- "loopless"
----- "regular"
----- "simple"
----- "slcRegular"
----- "ternary"
--------------------------------------------------------------------
-------------------------------------------------------------------- 
diffMatrix(ZZ,ZZ,String) := (n,r,s) -> (
    diffMatrix(readBasisFile(n,r,s), readBasisFile(n-1,r,s))
);


------------------- createDiffMatrixFile ---------------------------
--------------------------------------------------------------------
----- INPUT: (n,r,string)
-----
----- OUTPUT: File
-----
----- DESCRIPTION: Given integers n, r and a string sthis will attempt
----- to construct the matrix representing the differential
----- s_{n,r} ---> s_{n-1,r} where s is the subcomplex reffered to by
----- property s with n elements and rank r from scratch. This should only
----- be needed once and rarely needed for users.
-----
----- Note s must be one of the strings below:
----- "basis"
----- "binary"
----- "cographic"
----- "graphic"
----- "loopless"
----- "regular"
----- "simple"
----- "slcRegular"
----- "ternary"
--------------------------------------------------------------------
-------------------------------------------------------------------- 
createDiffMatrixFile = method();
createDiffMatrixFile (ZZ, ZZ, String) := (n, r, s) -> (
    basisList := rankedBasis(n,r);
    fileName := "diffs/"|toString(s)|"/diff_"|toString(s)|"_"|toString(n)|"_"|toString(r)|".txt";
    f := openOut fileName;
    f << toExternalString diffMatrix(n, r, s);
    f << endl;
    close(f)
    )


--------------------------- homologyRank ---------------------------
--------------------------------------------------------------------
----- INPUT: (n,r,string)
-----
----- OUTPUT: Number
-----
----- DESCRIPTION: Given integers n, r, and string s this will
----- comput the homology at
----- s_{n+1,r} ----> s_{n,r} ---> s_{n-1,r}
----- where s is the subcomplex referened by the property s.
-----
----- Note s must be one of the strings below:
----- "basis"
----- "binary"
----- "cographic"
----- "graphic"
----- "loopless"
----- "regular"
----- "simple"
----- "slcRegular"
----- "ternary"
--------------------------------------------------------------------
-------------------------------------------------------------------- 
homologyRank = method();
homologyRank(ZZ,ZZ,String) := (n,r,s) -> (
    Cnr := readBasisFile(n,r,s);
    --
    dOut := diffMatrix(Cnr, readBasisFile(n-1,r,s));
    dIn  := diffMatrix(readBasisFile(n+1,r,s), Cnr);
    --
    #Cnr - rank dOut - rank dIn
);



----------------------------- createBasisFile ----------------------
--------------------------------------------------------------------
----- INPUT: (n,r)
-----
----- OUTPUT: file
-----
----- DESCRIPTION: Given integers n, r this will attempt to
----- compute the basis for basis for the full matroid complex
----- on n elements and rank r from scratch. This should only
----- be needed once and rarely needed for users.
--------------------------------------------------------------------
-------------------------------------------------------------------- 
createBasisFile = method();
createBasisFile (ZZ, ZZ) := (n, r) -> (
    basisList := rankedBasis(n,r);
    fileName := "bases/basis/basis_"|toString(n)|"_"|toString(r)|".txt";
    f := openOut fileName;
    for M in basisList do (
	B := bases(M);
	E := toList groundSet(M);
	f << toExternalString {E, B};
	f << endl;
	);
    close(f)
    )

-------------------------- readBasisFile ---------------------------
--------------------------------------------------------------------
----- INPUT: (n,r,string)
-----
----- OUTPUT: list
-----
----- DESCRIPTION: Given integers n, r and a string sthis will attempt
----- to read in the basis for s_{n,r} where this is the subcomplex
----- with property given by ths string s on n elements with rank r
----- Note s must be one of the strings below:
----- "basis"
----- "binary"
----- "cographic"
----- "graphic"
----- "loopless"
----- "regular"
----- "simple"
----- "slcRegular"
----- "ternary"
--------------------------------------------------------------------
-------------------------------------------------------------------- 
readBasisFile = method();
readBasisFile (ZZ, ZZ, String) := (n, r, s) -> ( 
    fileName := "bases/"|toString(s)|"/"|toString(s)|"_"|toString(n)|"_"|toString(r)|".txt";
    if not fileExists fileName then {}
    else (
	fileLines := lines get fileName;
	apply(fileLines, l -> (
	    L := value l;
	    E := L#0;
	    B := L#1;
	    matroid(E,B)
	    ))
    )
    )

-------------------------- createSubBasisFile ----------------------
--------------------------------------------------------------------
----- INPUT: (n,r,F,s)
-----
----- OUTPUT: file
-----
----- DESCRIPTION: Given integers n, r a function F, and a strings s
----- this will attempt to create a basis for the subcomplex with
----- property on n elements and rank r from scratch by pruning
----- the full matroid complex on n elements and rank r.
-----
----- The function F must be a function that takes in a matroid M and
----- returns true if M has property s and false if M does not have s.
----- This should only be needed once and rarely needed for users.
--------------------------------------------------------------------
-------------------------------------------------------------------- 
createSubBasisFile = method();
createSubBasisFile (ZZ,ZZ,Thing, String) := (n, r, F, s) -> (
    bigBasis := readBasisFile(n,r,"basis");
    subBasis := delete(null,apply(bigBasis, M ->(
		if F(M) == true then M
	    )));
    fileName := "bases/"|toString(s)|"/"|toString(s)|"_"|toString(n)|"_"|toString(r)|".txt";
    f := openOut fileName;
    for M in subBasis do (
	B := bases(M);
	E := toList groundSet(M);
	f << toExternalString {E, B};
	f << endl;
	);
    close(f)
    )


----------------------------- isGraphic ---------------------------
--------------------------------------------------------------------
----- INPUT: (Matroid)
-----
----- OUTPUT: Boolean
-----
----- DESCRIPTION: Given a matroid checks whether it is graphic.
--------------------------------------------------------------------
-------------------------------------------------------------------- 
isGraphic = method(TypicalValue => Boolean);
isGraphic (Matroid) := (M) -> (
    -- Uses Tutte's minor characterization.
    U := uniformMatroid(2, 4);
    F7 := specificMatroid "fano";
    F7star := dual F7;
    K5 := completeGraph 5;
    K33 := completeMultipartiteGraph {3,3};
    MK5star := dual matroid K5;
    MK33star := dual matroid K33;
    not any({U, F7, F7star, MK5star, MK33star}, N -> hasMinor(M, N))
    )

----------------------------- isCographic --------------------------
--------------------------------------------------------------------
----- INPUT: (Matroid)
-----
----- OUTPUT: Boolean
-----
----- DESCRIPTION: Given a matroid checks whether it is cographic.
--------------------------------------------------------------------
-------------------------------------------------------------------- 
isCographic = method();
isCographic (Matroid) := (M) -> (
    isGraphic(dual(M))
    )

---------------------------- isRegularMatroid ----------------------
--------------------------------------------------------------------
----- INPUT: (Matroid)
-----
----- OUTPUT: Boolean
-----
----- DESCRIPTION: Given a matroid checks whether it is regular.
--------------------------------------------------------------------
-------------------------------------------------------------------- 
isRegularMatroid = method(TypicalValue => Boolean);
isRegularMatroid (Matroid) := (M) -> (
    U := uniformMatroid(2, 4);
    F7 := specificMatroid "fano";
    F7star := dual F7;
    not any({U, F7, F7star}, N -> hasMinor(M, N))
    )


---------------------------- isTernary -----------------------------
--------------------------------------------------------------------
----- INPUT: (Matroid)
-----
----- OUTPUT: Boolean
-----
----- DESCRIPTION: Given a matroid checks whether it is ternary.
--------------------------------------------------------------------
-------------------------------------------------------------------- 
isTernary = method(TypicalValue => Boolean);
isTernary (Matroid) := (M) -> (
    U25 := uniformMatroid(2, 5);
    U35 := uniformMatroid(3, 5);
    F7 := specificMatroid "fano";
    F7star := dual F7;
    not any({U25, U35, F7, F7star}, N -> hasMinor(M, N))
    )
--------------------------------------------------------------------
--------------------------------------------------------------------
----- Begining of the tests and the documentation
--------------------------------------------------------------------
--------------------------------------------------------------------

load ("./MatroidComplexes/tests.m2")
beginDocumentation()
load ("./MatroidComplexes/doc.m2")

end


--------------------------------------------------------------------
--------------------------------------------------------------------
----- Begining of sandbox
--------------------------------------------------------------------
--------------------------------------------------------------------

---
---
restart
uninstallPackage "MatroidComplexes"
restart
installPackage "MatroidComplexes"
check "MatroidComplexes"
installPackage "MatroidComplexes"



--------------------------------------------------------------------
--------------------------------------------------------------------
----- These are the runs to generate the basis files, display the
----- basis data as a table, and compute the differentials files
--------------------------------------------------------------------
--------------------------------------------------------------------

--------------------------  FULL BASIS   ---------------------------
--------------------------------------------------------------------
apply(toList(1..8), n->apply(toList(0..n), r->( print toString(n,r); time createBasisFile(n,r))));
apply({0,1,2,3,6,7,8,9}, r -> ( print toString(9,r); time createBasisFile(9,r)))
apply({4,5}, r -> ( print toString(9,r); time createBasisFile(9,r)))
--
basisFull = apply(toList(1..8), n->apply(toList(0..n), r->(readBasisFile(n,r,"basis"))))
apply(basisFull,n->apply(n,r->#r))
--
apply(toList(1..8), n->apply(toList(0..n), r->( print toString(n,r); time createDiffMatrixFile(n,r, "basis"))));
apply({0,1,2,3,6,7,8,9}, r -> ( print toString(9,r); time createDiffMatrixFile(9,r, "basis")));


-----------------------  GRAPHIC BASIS   ---------------------------
--------------------------------------------------------------------
apply(toList(1..8), n->apply(toList(0..n), r->( print toString(n,r); time createSubBasisFile(n,r,isGraphic,"graphic"))));
apply({0,1,2,3,6,7,8,9}, r -> ( print toString(9,r); time createSubBasisFile(9,r,isGraphic,"graphic")))
--
graphic = apply(toList(1..8), n->apply(toList(0..n), r->(readBasisFile(n,r,"graphic"))))
apply(graphic,n->apply(n,r->#r))
--
apply(toList(1..8), n->apply(toList(0..n), r->( print toString(n,r); time createDiffMatrixFile(n,r, "graphic"))));
apply({0,1,2,3,6,7,8,9}, r -> ( print toString(9,r); time createDiffMatrixFile(9,r, "graphic")));


-----------------------  COGRAPHIC BASIS   -------------------------
--------------------------------------------------------------------
apply(toList(1..8), n->apply(toList(0..n), r->( print toString(n,r); time createSubBasisFile(n,r,isCographic,"cographic"))));
apply({0,1,2,3,6,7,8,9}, r -> ( print toString(9,r); time createSubBasisFile(9,r,isCographic,"cographic")))
--
cographic = apply(toList(1..8), n->apply(toList(0..n), r->(readBasisFile(n,r,"cographic"))))
apply(cographic,n->apply(n,r->#r))
--
apply(toList(1..8), n->apply(toList(0..n), r->( print toString(n,r); time createDiffMatrixFile(n,r, "cographic"))));
apply({0,1,2,3,6,7,8,9}, r -> ( print toString(9,r); time createDiffMatrixFile(9,r, "cographic")));

-----------------------  REGULAR BASIS   -------------------------
--------------------------------------------------------------------
apply(toList(1..8), n->apply(toList(0..n), r->( print toString(n,r); time createSubBasisFile(n,r,isRegularMatroid,"regular"))));
apply({0,1,2,3,6,7,8,9}, r -> ( print toString(9,r); time createSubBasisFile(9,r,isRegularMatroid,"regular")))
--
regMatroids = apply(toList(1..8), n->apply(toList(0..n), r->(readBasisFile(n,r,"regular"))))
apply(regMatroids,n->apply(n,r->#r))
--
apply(toList(1..8), n->apply(toList(0..n), r->( print toString(n,r); time createDiffMatrixFile(n,r, "regular"))));
apply({0,1,2,3,6,7,8,9}, r -> ( print toString(9,r); time createDiffMatrixFile(9,r, "regular")));

-------------------------  BINARY BASIS  ---------------------------
--------------------------------------------------------------------
apply(toList(1..8), n->apply(toList(0..n), r->( print toString(n,r); time createSubBasisFile(n,r,isBinary,"binary"))));
apply({0,1,2,3,6,7,8,9}, r -> ( print toString(9,r); time createSubBasisFile(9,r,isBinary,"binary")))
--
binary = apply(toList(1..8), n->apply(toList(0..n), r->(readBasisFile(n,r,"binary"))))
apply(binary,n->apply(n,r->#r))
--
apply(toList(1..8), n->apply(toList(0..n), r->( print toString(n,r); time createDiffMatrixFile(n,r, "binary"))));
apply({0,1,2,3,6,7,8,9}, r -> ( print toString(9,r); time createDiffMatrixFile(9,r, "binary")));

------------------------  TERNARY BASIS  ---------------------------
--------------------------------------------------------------------
apply(toList(1..8), n->apply(toList(0..n), r->( print toString(n,r); time createSubBasisFile(n,r,isTernary,"ternary"))));
apply({0,1,2,3,6,7,8,9}, r -> ( print toString(9,r); time createSubBasisFile(9,r,isTernary,"ternary")))
--
ternary = apply(toList(1..8), n->apply(toList(0..n), r->(readBasisFile(n,r,"ternary"))))
apply(ternary,n->apply(n,r->#r))
--
apply(toList(1..8), n->apply(toList(0..n), r->( print toString(n,r); time createDiffMatrixFile(n,r, "ternary"))));
apply({0,1,2,3,6,7,8,9}, r -> ( print toString(9,r); time createDiffMatrixFile(9,r, "ternary")));

-------------------------  SIMPLE BASIS  ---------------------------
--------------------------------------------------------------------
apply(toList(1..8), n->apply(toList(0..n), r->( print toString(n,r); time createSubBasisFile(n,r,isSimple,"simple"))));
apply({0,1,2,3,6,7,8,9}, r -> ( print toString(9,r); time createSubBasisFile(9,r,isSimple,"simple")))
--
simple = apply(toList(1..8), n->apply(toList(0..n), r->(readBasisFile(n,r,"simple"))))
apply(simple,n->apply(n,r->#r))
--
apply(toList(1..8), n->apply(toList(0..n), r->( print toString(n,r); time createDiffMatrixFile(n,r, "simple"))));
apply({0,1,2,3,6,7,8,9}, r -> ( print toString(9,r); time createDiffMatrixFile(9,r, "simple")));


-------------------------  slcREGULAR BASIS  ---------------------------
--------------------------------------------------------------------
slcReg = apply(toList(1..15), n->apply(toList(0..n), r->(try readBasisFile(n,r,"slcRegular") else infinity)))
apply(slcReg,n->apply(n,r->(if not instance(r,List) then infinity else #r)))
--
apply(toList(1..15), n->apply(toList(0..n), r->( print toString(n,r); time createDiffMatrixFile(n,r, "slcRegular"))));

