dune-navier-stokes-two-phase
============================

[dune-navier-stokes-two-phase][0] is a
[Distributed and Unified Numerics Environment][1] module which implements a front
tracking finite element scheme to solve a single/two-phase incompressible
(Navier--)Stokes flow in 2D and 3D domains. It supports various finite
element spaces and boundary condition, several solvers and different
preconditioners. The full Navier--Stokes flow can be solved with a linearized
scheme or with a fixed-point iteration. It provides smoothing, remeshing and
ALE functionalities to overcome mesh distortions.

A detailed description of this method can be found in the paper

[M.~Agnese and R.~N\"urnberg, Fitted Finite Element Discretization of Two-Phase
{S}tokes Flow, Int.\ J.\ Numer.\ Meth.\ Fluids, 2016.][2].

**This is the paper I would ask everyone to cite when using
dune-navier-stokes-two-phase.**

License
-------

The dune-navier-stokes-two-phase library, headers and test programs are free
open-source software, licensed under version 2 or later of the GNU General
Public License.

See the file [LICENSE.md][3] for full copying permissions.

Installation
------------

For installation instructions please see the [DUNE website][4].

[0]: https://github.com/magnese/dune-repo/blob/master/dune-navier-stokes-two-phase/
[1]: https://www.dune-project.org/
[2]: http://onlinelibrary.wiley.com/doi/10.1002/fld.4237/abstract/
[3]: https://github.com/magnese/dune-repo/blob/master/dune-navier-stokes-two-phase/LICENSE.md
[4]: https://www.dune-project.org/doc/installation/
