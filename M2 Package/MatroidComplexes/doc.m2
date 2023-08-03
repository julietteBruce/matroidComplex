doc ///
   Key 
      MatroidComplexes
   Headline 
      Tools for computing the matroid chain complex
   Description
    Text
      The authors of package used a combination of high throughput
      high perfomance computing and sparse linear algebra to compute 
      the syzygies of $\mathbb P^1 x \mathbb P^1$ under the $D$-uple Veronese embedding
      for a number of values of $D=\{d_{1},d_{2}\}$. See the paper ``ADD AT SOME POINT'' by Bruce, 
      Corey, Erman, Goldstein, Laudone, Pirnes, and Yang, which we refer to
      as [BCEGLPY] throughout the documentation for this package.
      The goal of this package is to make this data more accessible and easy to 
      use by providing a way to access it via Macaulay2.

      Most functions have been implemented with three parameters $(B,D)$ where
      the goal is to compute the syzygies of the pushforward of
      the line bundle $\mathcal O_{\mathbb P^1 x \mathbb P^1}(B)$ under the $D$-uple
      embedding. Our hope is that as we (or others) are able to compute new data, we 
      will be able to update the package.

      
      One of the main functions is totalBettiTally which produces the standard
      Betti tables.  Other main functions refine the data in the Betti table
      by providing the multigraded Betti number or the Schur functor decompositions,
      or by computing statistics related to the Betti table (e.g. the BoijSoederberg
      coefficients) or related to the SchurFunctor decomposition.      

///


doc ///
   Key 
    withoutOddAut
    (withoutOddAut,Matroid)
   Headline
    determines if a matroid admits odd automorphisms
   Usage
    withoutOddAut(M)
   Inputs
    M: Matroid
   Outputs
    : Boolean
   Description
    Text
      Internally, a matroid $M$ on $n$ ground set elements is identified on the ground set 
      $\{0,1,2,\dots,n-1\}.$ Thus, any automorphism of $M$ is given by a permutation of 
      $\{0,1,2,\dots,n-1\}.$ Given a matroid $M$, this method checks whether the sign of any
      such permutation is negative, corresponding to an odd automorphism, and returns "true" 
      if a negative permutation is found, otherwise returns "false". 
      
      Those matroids which admit odd automorphisms vanish in the matroid complex. In the 
      example below, we ask whether the uniform matroids $U_{2,4}$ admits odd automorphisms. 
      That is, a matroid automorphism given by an odd permutation of its ground set $\{0,1,2,3\}.$
    Example
      M = uniformMatroid(2,4);
      withoutOddAut(M)

///

doc ///
   Key 
    rankedBasis
    (rankedBasis,Number,Number)
   Headline
    a list containing the basis elements of $C^r_n$
   Usage
    rankedBasis(n,r)
   Inputs
    n: Number
    r: Number
   Outputs
    : List
   Description
    Text
      Given a pair $(n,r)$, $0\leq r\leq n\leq n,$ this method returns a list of rank $r$ 
      matroids on $n$ ground set elements which do not admit odd automorphisms. These matroids
      constitute the basis of $C^r_n$ in the matroid complex.

      In example below, we generate the basis of rank $3$ matroids on $6$ ground set elements. 
    Example
      rankedBasis(6,3)
     
///

doc ///
   Key 
    diffMatrixColumn
    (diffMatrixColumn,Matroid,List)
   Headline
    the column vector which is the image of a basis under the deletion differential with 
    respect to standard bases
   Usage
    diffMatrixColumn(M,T)
   Inputs
    M: Matroid
    T: List
   Outputs
    : List
   Description
    Text
      Given a matroid representing a basis element of (a subspace of) C_n, return the corresponding column vector of the
      matrix representing the deletion differential from C_n to C_{n-1}
      with respect to the standard bases on Q^{dim C_n} and
      Q^{dim C_{n-1}}.
      
      In example below, we begin by defining the cyclic matroid $M$ associated to the wheel graph $W_3$ (with three spokes), 
      which is a element of the basis of $C^3_6,$ as can be verified by generating rankedBasis(6,3). Then, we generate the 
      column vector associated to $M$ under the deletion differential from $C_6$ to $C_5$, with respect to the standard bases 
      on $Q^{dim C_n}$ and $Q^{dim C_{n-1}}$.
    Example
      M = wheel 3;
      targetBasis = rankedBasis(6,3);
      diffMatrixColumn(M,targetBasis)
  
///

doc ///
   Key 
    diffMatrix
    (diffMatrix,List,List)
   Headline
    the matrix representing the deletion differential with respect to standard bases
   Usage
    diffMatrix(S,T)
   Inputs
    S: List
    T: List
   Outputs
    : Matrix
   Description
    Text
      Given the pair of bases for a subspace $V$ of $C_n$ and a basis for a subspace $W$
      of $C_{n-1}$ containing the image of $V$ under the deletion differential, return the 
      matrix representing the standard bases on $QQ^{dim V}$ and $QQ^{dim W}$, respectively. 
      The first element of the first list will correspond to the basis vector 
      $e_1 = (1,0,...,0)^T$, etc.
 
      In the following example, we generate the matrix associated to the deletion differential 
      from $C^3_7$ to $C^3_6$ with respect to the standard bases $QQ^9$ and $QQ^2,$ the dimensions
      of which we verified by running rankedBasis(7,3) and rankedBasis(6,3), respectively. 
    Example
      sourceBasis = rankedBasis(7,3);
      targetBasis = rankedBasis(6,3);
      diffMatrix(sourceBasis,targetBasis)
///
