hemco
=====

This is the Harvard-NASA Emissions Component (HEMCO) source code. The Harvard-NASA Emissions Component (HEMCO) is a software
component for computing (atmospheric) emissions from different sources, regions, and species on a user-defined grid. It
can combine, overlay, and update a set of data inventories ('base emissions') and scale factors, as specified by the
user through the HEMCO configuration file. Emissions that depend on environmental variables and non-linear
parameterizations are calculated in separate HEMCO extensions. HEMCO can be run in standalone mode or coupled to an
atmospheric model.
A more detailed description of HEMCO is given in Keller et al. (2014).

This is just the plain HEMCO code without an explicit model interface. For the standalone HEMCO model code, see
github.com/christophkeller/hemco_standalone_full.

References: C. A. Keller, M. S. Long, R. M. Yantosca, A. M. Da Silva, S. Pawson, D. J. Jacob: HEMCO v1.0: a versatile,
ESMF-compliant component for calculation emissions in atmospheric models. Geosci. Model Dev., 7, 1409-1417, 2014.

Contact: ckeller@seas.harvard.edu
