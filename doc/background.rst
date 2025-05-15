Background
==========

History
-------

The development of AMICI started out in 2015 by Fabian Fröhlich in the
Jan Hasenauer lab as an *Advanced MATLAB Interface to CVODES and IDAS*.
The goal was to provide a fast and flexible interface to the SUNDIALS
:term:`CVODES` and :term:`IDAS` solvers for the efficient simulation and
sensitivity analysis of large-scale biochemical reaction networks while
supporting high-level model descriptions in :term:`SBML` as the de facto
standard for systems biology models.
In contrast to the (no longer maintained)
`sundialsTB <https://computing.llnl.gov/projects/sundials/sundials-software>`__
MATLAB interface, all necessary functions were transformed into native
C++ code, which allowed for significantly faster simulations.

Over time, various other researchers contributed to AMICI's development.
AMICI's C core was reshaped into a reusable stand-alone C++ library,
that was the basis for its Python interface. These additional interfaces lead
to the reinterpretation of "AMICI" as "Advanced *Multi-language* Interface to
CVODES and IDAS".

Due to changes in MATLAB's symbolic math library and to foster open and
reusable research, the development of AMICI's MATLAB interface
was meanwhile discontinued, and the focus shifted to the Python interface.

Over the years, AMICI has evolved into a powerful and flexible
tool for the simulation and sensitivity analysis of large-scale
biochemical reaction networks, supporting different input formats,
different types of sensitivity analysis, discrete events,
and many other features.
Meanwhile, the use of AMICI has been documented in around 100 research papers
from different labs and companies.

AMICI is currently maintained by Fabian Fröhlich
(`@FFroehlich <https://github.com/FFroehlich>`__,
`Fröhlich lab <https://www.frohlichlab.com/>`__) and Daniel Weindl
(`@dweindl <https://github.com/dweindl>`__,
`Hasenauer lab <https://www.mathematics-and-life-sciences.uni-bonn.de/en/research/hasenauer-group>`__).

Publications on various features of AMICI
-----------------------------------------

Some mathematical background for AMICI is provided in the following
publications:

* Fröhlich, F., Kaltenbacher, B., Theis, F. J., & Hasenauer, J. (2017).
  **Scalable Parameter Estimation for Genome-Scale Biochemical Reaction Networks.**
  PLOS Computational Biology, 13(1), e1005331.
  doi:`10.1371/journal.pcbi.1005331 <https://doi.org/10.1371/journal.pcbi.1005331>`_.

* Fröhlich, F., Theis, F. J., Rädler, J. O., & Hasenauer, J. (2017).
  **Parameter estimation for dynamical systems with discrete events and logical
  operations.** Bioinformatics, 33(7), 1049-1056.
  doi:`10.1093/bioinformatics/btw764 <https://doi.org/10.1093/bioinformatics/btw764>`_.

* Terje Lines, Glenn, Łukasz Paszkowski, Leonard Schmiester, Daniel Weindl,
  Paul Stapor, and Jan Hasenauer. 2019. **Efficient Computation of Steady States
  in Large-Scale ODE Models of Biochemical Reaction Networks.**
  *IFAC-PapersOnLine* 52 (26): 32–37.
  DOI: `10.1016/j.ifacol.2019.12.232 <https://doi.org/10.1016/j.ifacol.2019.12.232>`_.

* Stapor, Paul, Fabian Fröhlich, and Jan Hasenauer. 2018.
  **Optimization and Profile Calculation of ODE Models Using Second Order
  Adjoint Sensitivity Analysis.** *Bioinformatics* 34 (13): i151–i159.
  DOI: `10.1093/bioinformatics/bty230 <https://doi.org/10.1093/bioinformatics/bty230>`_.

* Lakrisenko, Polina, Paul Stapor, Stephan Grein, Łukasz Paszkowski,
  Dilan Pathirana, Fabian Fröhlich, Glenn Terje Lines, Daniel Weindl,
  and Jan Hasenauer. 2023.
  **Efficient Computation of Adjoint Sensitivities at Steady-State in ODE Models
  of Biochemical Reaction Networks.** *PLoS Comput Biol* 19(1): e1010783.
  DOI: `10.1371/journal.pcbi.1010783 <https://doi.org/10.1371/journal.pcbi.1010783>`_.

* L. Contento, P. Stapor, D. Weindl, and J. Hasenauer. 2023.
  **A more expressive spline representation for SBML models improves code generation performance in AMICI**,
  In: Pang, J., Niehren, J. (eds) Computational Methods in Systems Biology.
  CMSB 2023. *Lecture Notes in Computer Science*, vol 14137. Springer, Cham.
  DOI: `10.1007/978-3-031-42697-1_3 <https://doi.org/10.1007/978-3-031-42697-1_3>`_.
  Preprint available at `bioRxiv <https://doi.org/10.1101/2023.06.29.547120>`_.

* Lakrisenko, Polina, Dilan Pathirana, Daniel Weindl, and Jan Hasenauer. 2024.
  **Exploration of methods for computing sensitivities in ODE models at dynamic and steady states.** *arXiv:2405.16524 [q-bio.QM]*.
  DOI: `10.48550/arXiv.2405.16524 <https://doi.org/10.48550/arXiv.2405.16524>`_.


.. note::

   Implementation details of the latest AMICI versions may differ from the ones
   given in the referenced manuscripts.


Third-Party numerical algorithms used by AMICI
----------------------------------------------

AMICI uses the following packages from `SUNDIALS <https://computing.llnl.gov/projects/sundials>`__:

* CVODES:

  The sensitivity-enabled ODE solver in SUNDIALS. Radu Serban
  and Alan C. Hindmarsh. *ASME 2005 International Design Engineering
  Technical Conferences and Computers and Information in Engineering
  Conference*. American Society of Mechanical Engineers, 2005.
  `PDF <http://proceedings.asmedigitalcollection.asme.org/proceeding.aspx?articleid=1588657>`__

* IDAS

AMICI uses the following packages from `SuiteSparse <https://github.com/DrTimothyAldenDavis/SuiteSparse>`__:

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

Others:

* SuperLU_MT

  "A general purpose library for the direct solution of large,
  sparse, nonsymmetric systems of linear equations"
  (https://crd-legacy.lbl.gov/~xiaoye/SuperLU/#superlu_mt).
  SuperLU_MT is optional and is so far only available from the C++ interface.
