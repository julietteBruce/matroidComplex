
--------------------------------------------------------------------
--------------------------------------------------------------------
----- Tests for withoutOddAut.
--------------------------------------------------------------------
--------------------------------------------------------------------
TEST ///
    M = uniformMatroid(2,4);
    assert (withoutOddAut(M) == false)
///

TEST ///
    M = uniformMatroid(0,4);
    assert (withoutOddAut(M) == false)
///

TEST /// 
    M = wheel 3;
    assert (withoutOddAut(M) == true)
///

--------------------------------------------------------------------
--------------------------------------------------------------------
----- Tests for rankedBasis.
--------------------------------------------------------------------
--------------------------------------------------------------------
TEST ///
    B44 = rankedBasis(4,4);
    assert (B44 == {})
///

TEST ///
    B63 = rankedBasis(6,3);
    assert (B63 == {Matroid{bases ⇒ {set {0,1,3},set {0,2,3},set {1,2,3},set {0,1,4},set {0,2,4},set {1,2,4},set {1,3,4},set {2,3,4},set {0,1,5},set {0,2,5},set {1,2,5},set {0,3,5},set {2,3,5},set {0,4,5},set {1,4,5},set {2,4,5},set {3,4,5}},cache ⇒ (CacheTable{…3…}),groundSet ⇒ set {0,1,2,3,4,5},rank ⇒ 3},Matroid{bases ⇒ {set {0,1,3},set {0,2,3},set {1,2,3},set {0,1,4},set {0,2,4},set {1,2,4},set {1,3,4},set {2,3,4},set {0,1,5},set {0,2,5},set {1,2,5},set {0,3,5},set {2,3,5},set {0,4,5},set {1,4,5},set {3,4,5}},cache ⇒ (CacheTable{…3…}),groundSet ⇒ set {0,1,2,3,4,5},rank ⇒ 3}})
///

--------------------------------------------------------------------
--------------------------------------------------------------------
----- Tests for diffMatrixColumn.
--------------------------------------------------------------------
--------------------------------------------------------------------

TEST (///
    M = uniformMatroid(7,3);
    targetBasis = rankedBasis(6,3);
    C = diffMatrixColumn(M,targetBasis);
    assert (C == {0,0})
///

--------------------------------------------------------------------
--------------------------------------------------------------------
----- Tests for diffMatrix.
--------------------------------------------------------------------
--------------------------------------------------------------------
TEST ///
    sourceBasis = rankedBasis(7,3);
    targetBasis = rankedBasis(6,3);
    M = diffMatrix(sourceBasis,targetBasis);
    assert (M == matrix{{-1,-1,0,0,0,0,-1,-1,0},{0,0,-1,-1,-1,-1,1,0,-1}})

///




