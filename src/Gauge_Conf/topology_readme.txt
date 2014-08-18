In the continuum the topological charge is defined by

Q=\frac{1}{32 \pi^2} \int \epsilon_{\mu\nu\rho\sigma}Tr(F_{\mu\nu}F_{\rho\sigma}) d^x 

The discretization used on the lattice is obtained by using the clover form of the 
discretized field-strenth, i.e. the clover Q_{\mu\nu} ("quadrifoglio" in the code) is 
given by

Q_{\mu\nu} = \exp(4 i a^2 F_{\nu\nu ) ~ 1 + 4 i a^2 F_{\mu\nu}

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
