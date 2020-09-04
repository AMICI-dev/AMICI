Background
==========

*This section is to be extended.*

Publications on various features of AMICI
-----------------------------------------

Some mathematical background for AMICI is provided in the following
publications:

* Fröhlich, F., Kaltenbacher, B., Theis, F. J., & Hasenauer, J. (2017).
  Scalable Parameter Estimation for Genome-Scale Biochemical Reaction Networks.
  PLOS Computational Biology, 13(1), e1005331.
  doi:`10.1371/journal.pcbi.1005331 <https://doi.org/10.1371/journal.pcbi.1005331>`_.

* Fröhlich, F., Theis, F. J., Rädler, J. O., & Hasenauer, J. (2017).
  Parameter estimation for dynamical systems with discrete events and logical
  operations. Bioinformatics, 33(7), 1049-1056.
  doi:`10.1093/bioinformatics/btw764 <https://doi.org/10.1093/bioinformatics/btw764>`_.

* Terje Lines, Glenn, Łukasz Paszkowski, Leonard Schmiester, Daniel Weindl,
  Paul Stapor, and Jan Hasenauer. 2019. "Efficient Computation of Steady States
  in Large-Scale Ode Models of Biochemical Reaction Networks.
  *IFAC-PapersOnLine* 52 (26): 32–37.
  DOI: `10.1016/j.ifacol.2019.12.232 <https://doi.org/10.1016/j.ifacol.2019.12.232>`_.

* Stapor, Paul, Fabian Fröhlich, and Jan Hasenauer. 2018.
  "Optimization and Profile Calculation of ODE Models Using Second Order
  Adjoint Sensitivity Analysis." *Bioinformatics* 34 (13): i151–i159.
  DOI: `10.1093/bioinformatics/bty230 <https://doi.org/10.1093/bioinformatics/bty230>`_.


.. note::

   Implementation details of the latest AMICI versions may differ from the ones
   given in the references manuscripts.


Third-Party numerical algorithms used by AMICI
----------------------------------------------

AMICI uses the following packages from SUNDIALS:

* CVODES:

  The sensitivity-enabled ODE solver in SUNDIALS. Radu Serban
  and Alan C. Hindmarsh. *ASME 2005 International Design Engineering
  Technical Conferences and Computers and Information in Engineering
  Conference*. American Society of Mechanical Engineers, 2005.
  `PDF <http://proceedings.asmedigitalcollection.asme.org/proceeding.aspx?articleid=1588657>`__

* IDAS

AMICI uses the following packages from SuiteSparse:

* Algorithm 907: **KLU** A Direct Sparse Solver for Circuit Simulation
  Problems. Timothy A. Davis, Ekanathan Palamadai Natarajan,
  *ACM Transactions on Mathematical Software*, Vol 37, Issue 6, 2010,
  pp 36:1-36:17. `PDF <http://dl.acm.org/authorize?305534>`__

* Algorithm 837: **AMD**, an approximate minimum degree ordering
  algorithm, Patrick R. Amestoy, Timothy A. Davis, Iain S. Duff,
  *ACM Transactions on Mathematical Software*, Vol 30, Issue 3, 2004,
  pp 381-388. `PDF <http://dl.acm.org/authorize?733169>`__

* Algorithm 836: **COLAMD**, a column approximate minimum degree ordering
  algorithm, Timothy A. Davis, John R. Gilbert, Stefan I. Larimore,
  Esmond G. Ng *ACM Transactions on Mathematical Software*, Vol 30,
  Issue 3, 2004, pp 377-380. `PDF <http://dl.acm.org/authorize?734450>`__
