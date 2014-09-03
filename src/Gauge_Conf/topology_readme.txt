Topology for the general group: in the continuum the topological charge is defined by

Q=\frac{1}{64 \pi^2} \int \epsilon_{\mu\nu\rho\sigma} F_{\mu\nu}^a F_{\rho\sigma}^a d^x

with the normalization Tr(T_iT_j)=K \delta_{ij} such that the longest root of the
representation is normalized to 1.
[see C. W. Bernard, N. H. Christ, A. H. Guth, E. J. Weinberg Phys. Rev. D 16, 2967 (1977) ]

For SU(N) this bases is that of the generalized Gell-Mann matrices divided by two 
and in this representation Tr(T^aT^b)=(1/2)\delta^{ab}, so that 

F_{\mu\nu}^aF_{\rho\sigma}^a=2Tr(F_{\mu\nu}F_{\rho\sigma})

and thus for SU(N)

Q=\frac{1}{32 \pi^2} \int \epsilon_{\mu\nu\rho\sigma} F_{\mu\nu}^a F_{\rho\sigma}^a d^x

For G2 a basis of the algebra is explicitly constructed in 
S. L. Cacciatori, B. L. Cerchiai, A. Della Vedova, G. Ortenzi, 
A. Scotti J. Math. Phys. 46, 083512 (2005)
and imposing as normalization that the longest root is equal to 1, we get 
Tr(T^aT^b)=\delta^{ab}, so G_2 we have

Q=\frac{1}{64 \pi^2} \int \epsilon_{\mu\nu\rho\sigma} Tr (F_{\mu\nu}F_{\rho\sigma}) d^x

-------------
| FOR SU(N) |
-------------

The discretization used on the lattice is obtained by using the clover form of the 
discretized field-strenth, i.e. the clover Q_{\mu\nu} ("quadrifoglio" in the code) is 
given by

Q_{\mu\nu}  ~ 4 + 4 i a^2 F_{\mu\nu}

and thus

F_{\mu\nu} = \frac{1}{8i}( Q_{\mu\nu} - Q_{\mu\nu}^{dag} )  

By using this expression we get

-Tr(F_{\mu\nu}F_{\rho\sigma}) = 
          =\frac{1}{2^5} ReTr(Q_{\mu\nu}[Q_{\rho\sigma}-Q_{\rho\sigma}^{dag}])

We now note that of the 24 terms \epsilon_{\mu\nu\rho\sigma}F_{\mu\nu}F_{\rho\sigma}
only 3 are independent:

2: single exchange on the first term   F_{\mu\nu}F_{\rho\sigma} -> F_{\nu\mu}F_{\rho\sigma}
2: single exchange on the second term  F_{\mu\nu}F_{\rho\sigma} -> F_{\mu\nu}F_{\sigma\rho}
2: exchage                             F_{\mu\nu}F_{\rho\sigma} -> F_{\rho\sigma}F_{\mu\nu} 
3: cyclical permutations of indices

thus 

\frac{1}{32} \epsilon_{\mu\nu\rho\sigma}Tr(F_{\mu\nu}F_{\rho\sigma})=

  =-\frac{1}{32} * \frac{1}{2^5} * 2^3 * [ sum on cyclic permutations
     of \epsilon_{\mu\nu\rho\sigma} ReTr(Q_{\mu\nu}[Q_{\rho\sigma}-Q_{\rho\sigma}^{dag}]) ]
  
  =-\frac{1}{128} [ sum on cyclic permutations
     of \epsilon_{\mu\nu\rho\sigma} ReTr(Q_{\mu\nu}[Q_{\rho\sigma}-Q_{\rho\sigma}^{dag}]) ]

This is the expression used in the code. 

----------
| FOR G2 |
---------

For G2, the final expression for SU(N) has to be divided by 2
