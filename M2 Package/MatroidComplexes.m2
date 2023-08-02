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
        {Name => "Juliette Bruce",           Email => "juliette.bruce@berkeley.edu",       HomePage => "https://juliettebruce.github.io"},
        {Name => "Benjamin Ashlock",           Email => "bak6t@missouri.edu"},
        {Name => "Jacob Bucciarelli",             Email => "jbucciarelli@ksu.edu"},	     
        {Name => "Bailee Zacovic",          Email => "bzacovic@umich.edu"},
	},
  DebuggingMode => true,
  AuxiliaryFiles => true
  )

export {
  "withoutOddAut", 
  "rankedBasis", 
  "diffMatrixColumn", 
  "diffMatrix", 

  }

--------------------------------------------------------------------
--------------------------------------------------------------------
----- CODE
--------------------------------------------------------------------
--------------------------------------------------------------------

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
rankedBasis(Number,Number) := (n,r) -> (
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
viewHelp SchurVeronese
