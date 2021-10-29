(sec:example:subduction:3d:organization)=
# Organization of Simulation Parameters

PyLith automatically reads `pylithapp.cfg` from the current directory, if it exists.
As a result, we generally put all parameters common to a set of examples in this file to avoid duplicating parameters across multiple files.
Because we often use a single mesh for multiple simulations in a directory, we place all parameters related to our mesh and identifying the materials in our mesh in `pylithapp.cfg`.
We assign the bulk constitutive model and its parameters to each material in other files, because we vary those across the simulations.
In general, we place roller boundary conditions (Dirichlet boundary conditions constraining the degrees of freedom perpendicular to the boundary) on the lateral and bottom boundaries, so we include those in `pylithapp.cfg`.
In some simulations we will overwrite the values for parameters will values specific to a given example.
This file is also a convenient place to put basic solver parameters and to turn on Pyre journals for displaying informational and debugging messages.

The settings contained in `pylithapp.cfg` include:

**pylithapp.journal.info** Settings that control the verbosity of the output written to stdout for the different components.

**pylithapp.mesh_generator** Parameters for the type of mesh importer (generator), reordering of the mesh, and the mesh coordinate system.

**pylithapp.problem.materials** Basic parameters for each of the four materials, including the label, block id in the mesh file, discretization, and output writer.

**pylithapp.problem.bc** Parameters for Dirichlet boundary conditions on the lateral and bottom boundaries of the domain.

**pylithapp.problem.formulation.output** Settings related output of the solution over the domain and subdomain (ground surface).

**pylithapp.petsc** PETSc solver and logging settings.

### Coordinate system

We generated the mesh in a Cartesian coordinate system corresponding to a transverse Mercator projection.
We specify this geographic projection coordinate system in the `pylithapp.cfg` file, so that we can use other convenient georeferenced coordinate systems in the spatial databases.
PyLith will automatically transform points between compatible coordinate systems.
Our spatialdata library uses [Proj](https://proj.org) for geographic coordinate transformations, so we specify the projection using Proj syntax in the `crs_string` property:

```{code-block} cfg
---
caption: Excerpt from `pylithapp.cfg`
---
[pylithapp.mesh_generator.reader]
coordsys = spatialdata.geocoords.CSGeo
coordsys.space_dim = 3
coordsys.crs_string = +proj=tmerc +datum=WGS84 +lon_0=-122.6765 +lat_0=45.5231 +k=0.9996 +units=m
```

## Materials

The finite-element mesh marks cells for each material and associates an integer with each material. We specify this information in the `pylithapp.cfg` file and avoid duplicating it in each simulation parameter file.
To set up the materials, we first create an array of materials that defines the name for each material component.

```{code-block} cfg
---
caption: Excerpt from `pylithapp.cfg`
---
[pylithapp.problem]
materials = [slab, wedge, crust, mantle]

[pylithapp.problem.materials.slab]
label = Subducting slab
id = 1

[pylithapp.problem.materials.wedge]
label = Accretionary wedge
id = 2

[pylithapp.problem.materials.mantle]
label = Mantle
id = 3

[pylithapp.problem.materials.crust]
label = Continental crust
id = 4
```

In this set of examples, we will consider cases in which all materials are linear, isotropic elastic and cases where the crust and wedge are linear, isotropic elastic but the slab and mantle are linear Maxwell viscoelastic.
As a result, we put the parameters for these two cases in separate `cfg` files with `mat_elastic.cfg` for the case with purely elastic models and `mat_viscoelastic.cfg` for the case with a mix of elastic and viscoelastic models.
Each of these files specifies the bulk constitutive model and spatial database to use for the properties for each material.
The values for the material properties are loosely based on a 3-D seismic velocity model for the Pacific Northwest {cite}`Stephenson:2007`.

### Boundary Conditions

For the Dirichlet boundary conditions, we specify the degree of freedom constrained, the name of the nodeset in the ExodusII file from CUBIT/Trelis that defines the boundary, and a label for the spatial database (required for informative error messages).
These settings constrain the y-displacement on the north (+y) boundary:

```{code-block} cfg
---
caption: Excerpt from `pylithapp.cfg`
---
[pylithapp.problem.bc.bc_ypos]
constrained_dof = [1]
label = boundary_ypos
db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.label = Dirichlet BC on +y
```

## Solver Parameters

We group solver parameters into a few different files to handle different cases.
The `pylithapp.cfg` contains tolerance values for the linear and nonlinear solvers and turns on some simple diagnostic information.
The file also directs PyLith to use an algebriac multigrid preconditioner, which works well for most problms that do not include a fault.
It will not work for problems with a fault.

```{code-block} cfg
---
caption: Excerpt from `pylithapp.cfg`
---
[pylithapp.petsc]
malloc_dump = true

ts_type = beuler
pc_type = ml

# Convergence parameters.
ksp_rtol = 1.0e-10
ksp_atol = 1.0e-12
ksp_max_it = 1000
ksp_gmres_restart = 50
ksp_error_if_not_converged = true

snes_rtol = 1.0e-10
snes_atol = 1.0e-10
snes_max_it = 3
snes_error_if_not_converged = true

# Monitors for debugging
ts_monitor = true
ksp_monitor = true
ksp_converged_reason = true
snes_monitor = true
snes_converged_reason = true
snes_linesearch_monitor = true
ksp_view = true
snes_view = true

# PETSc summary -- useful for performance information.
log_view = true
```

For simulations with a fault `solver_fieldsplit.cfg` provides settings for applying the algebraic multigrid preconditioner to the elasticity portion of the system Jacobian matrix and a preconditioner suitable for our fault formulation to the Lagrange multiplier portion.
