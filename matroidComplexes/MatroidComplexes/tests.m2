
--------------------------------------------------------------------
--------------------------------------------------------------------
----- Tests for makeBettiTally.
--------------------------------------------------------------------
--------------------------------------------------------------------
TEST ///
    assert (0 == 0 )
///


--------------------------------------------------------------------
--------------------------------------------------------------------
----- Tests for multiBettiP1P1.
--------------------------------------------------------------------
--------------------------------------------------------------------
TEST ///
    H = multiBettiP1P1({0,0},{2,4});
    f = t_0^14*t_1^14*t_2^29*t_3^27+t_0^14*t_1^14*t_2^28*t_3^28+t_0^14*t_1^14*t_2^27*t_3^29;
    assert (sub(H#(12,2),ring f)  == f)
///


--------------------------------------------------------------------
--------------------------------------------------------------------
----- Tests for schurBettiP1P1.
--------------------------------------------------------------------
--------------------------------------------------------------------
TEST ///
    H = schurBettiP1P1({0,0},{2,4});
    assert (H#(11,2) == {({14, 12, 28, 24}, 1), ({14, 12, 27, 25}, 1)})
///



--------------------------------------------------------------------
--------------------------------------------------------------------
----- Tests for totalBettiP1P1.
--------------------------------------------------------------------
--------------------------------------------------------------------
TEST ///
    H = totalBettiP1P1({0,0},{2,4});
    assert (H#(1,1) == 75 )
///


--------------------------------------------------------------------
--------------------------------------------------------------------
----- Tests for totalBettiTallyP1P1.
--------------------------------------------------------------------
--------------------------------------------------------------------
TEST ///
    assert (0 == 0 )
///


--------------------------------------------------------------------
--------------------------------------------------------------------
----- Tests for repsWithoutMultiplicityP1P1.
--------------------------------------------------------------------
--------------------------------------------------------------------
TEST ///
    G = schurBettiP1P1({0,0},{2,4});
    H = repsWithoutMultiplicityP1P1 G;
    assert (H#(11,2) == {{14, 12, 28, 24}, {14, 12, 27, 25}})
///


--------------------------------------------------------------------
--------------------------------------------------------------------
----- Tests for numDistinctRepsBettiP1P1.
--------------------------------------------------------------------
--------------------------------------------------------------------
TEST ///
    H = numDistinctRepsBettiP1P1({0,0},{2,4});
    assert (H#(11,2) == 2 )
///


--------------------------------------------------------------------
--------------------------------------------------------------------
----- Tests for numRepsBettiP1P1.
--------------------------------------------------------------------
--------------------------------------------------------------------
TEST ///
    H = numRepsBettiP1P1({0,0},{2,4});
    assert (H#(11,2) == 2)
///
