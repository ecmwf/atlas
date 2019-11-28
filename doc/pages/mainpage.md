Atlas                             {#mainpage}
=====

* *A library for numerical weather prediction and climate modelling*

This is a documentation for the Atlas project. To see high-level feature overview, or project goals, head over to the [project homepage](https://confluence.ecmwf.int/display/ATLAS). For different versions of this documentation see [here](http://download.ecmwf.int/test-data/atlas/docs).

@m_class{m-note m-info}

@parblock
@m_class{m-noindent}

Atlas is a ECMWF library for parallel data-structures supporting unstructured
grids and function spaces, with the aim to investigate alternative more scalable
dynamical core options for Earth System models, and to support modern interpolation
and product generation software
@endparblock

Atlas is predominantly C++ code, with main features available to Fortran codes
through a F2003 interface. It requires some flavour of Unix (such as Linux).
It is known to run on a number of systems, some of which are directly supported
by ECMWF.


What's new?                  {#mainpage-whats-new}
-----------

Curious about what was added or improved recently? Check out the
[Changelog](https://github.com/ecmwf/atlas/blob/master/CHANGELOG.md)

Getting started             {#mainpage-getting-started}
---------------

The best way to get started is to read the @ref getting-started guide.

Contributing                {#mainpage-contributing}
------------

Contributions to Atlas are welcome. In order to do so, please open an issue
where a feature request or bug can be discussed. Then issue a pull request
with your contribution. Pull requests must be issued against the develop branch.

Citing Atlas                {#mainpage-citing-atlas}
------------

If you publish work which mentions Atlas, or Atlas has been useful in your research,
please cite the following paper:

```bibtex
@article{DECONINCK2017188,
title = "Atlas : A library for numerical weather prediction and climate modelling",
journal = "Computer Physics Communications",
volume = "220",
pages = "188 - 204",
year = "2017",
issn = "0010-4655",
doi = "https://doi.org/10.1016/j.cpc.2017.07.006",
url = "http://www.sciencedirect.com/science/article/pii/S0010465517302138",
author = "Willem Deconinck and Peter Bauer and Michail Diamantakis and Mats Hamrud and Christian Kühnlein and Pedro Maciel and Gianmarco Mengaldo and Tiago Quintino and Baudouin Raoult and Piotr K. Smolarkiewicz and Nils P. Wedi",
keywords = "Numerical weather prediction, Climate, Earth system, High performance computing, Meteorology, Flexible mesh data structure"
}
```

