(sec:example:subduction:3d:step01)=
# Step 1: Axial Compression

We start with a very simple example of axial compression in the east-west direction with purely elastic material properties, and no faults ({numref}`fig:example:subduction:3d:step01:diagram`).
We impose axial compression using Dirichlet boundary conditions on the east (+x) and west (-x) boundaries and confine the domain in the north-south direction via zero displacement Dirichlet boundary conditions on the north (+y) and south (-y) boundaries.
We constrain the vertical displacement by imposing zero displacement boundary conditions on the bottom (-z) boundary.

:::{figure-md} fig:example:subduction:3d:step01:diagram
 <img src="figs/subduction3d_step01_diagram.*" alt="Diagram of Step 1 - Axial compression. This static simulation uses Dirichlet boundary conditions with axial compression in the east-west (x-direction), roller boundary conditions on the north, south, and bottom boundaries, and purely elastic properties." width="100%"/>

Diagram of Step 1 - Axial compression. This static simulation uses Dirichlet boundary conditions with axial compression in the east-west (x-direction), roller boundary conditions on the north, south, and bottom boundaries, and purely elastic properties.
:::

The `pylithapp.cfg` file creates an array of five boundary conditions, which impose zero displacements by default.
We overwrite this behavior in the `step01_axialdisp.cfg` file for the -x and +x boundaries using spatial databases with a single uniform displacement value to create the axial compression:

```{code-block} cfg
---
caption: Excerpt from `step01_axialdisp.cfg` showing the Dirichlet boundary condition settings.
---
# -x face
[pylithapp.problem.bc.x_neg]
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.label = Dirichlet BC on -x
db_auxiliary_field.values = [initial_amplitude_x, initial_amplitude_y, initial_amplitude_z]
db_auxiliary_field.data = [+2.0*m, 0.0*m, 0.0*m]

# +x face
[pylithapp.problem.bc.x_pos]
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.label = Dirichlet BC on +x
db_auxiliary_field.values = [initial_amplitude_x, initial_amplitude_y, initial_amplitude_z]
db_auxiliary_field.data = [-2.0*m, 0.0*m, 0.0*m]
```

As discussed in {ref}`sec:example:subduction:3d:organization`, we use `mat_elastic.cfg` to specify the parameters associated with linear, isotropic elastic bulk constitutive models for all of the materials for convenient reuse across several different simulations.

```{code-block} cfg
---
caption: Excerpt from `mat_elastic.cfg`
---
[pylithapp.problem.materials]
slab.bulk_rheology = pylith.materials.IsotropicLinearElasticity
wedge.bulk_rheology = pylith.materials.IsotropicLinearElasticity
crust.bulk_rheology = pylith.materials.IsotropicLinearElasticity
mantle.bulk_rheology = pylith.materials.IsotropicLinearElasticity

# Slab
[pylithapp.problem.materials.slab]
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.label = Properties for subducting slab
db_auxiliary_field.iohandler.filename = spatialdb/mat_slab_elastic.spatialdb

observers.observer.data_fields = [displacement, cauchy_stress, cauchy_strain]
observers.observer.trigger.num_skip = 1

auxiliary_subfields.density.basis_order = 0
bulk_rheology.auxiliary_subfields.bulk_modulus.basis_order = 0
bulk_rheology.auxiliary_subfields.shear_modulus.basis_order = 0

derived_subfields.cauchy_strain.basis_order = 1
derived_subfields.cauchy_stress.basis_order = 1


# Wedge
[pylithapp.problem.materials.wedge]
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.label = Properties for accretionary wedge
db_auxiliary_field.iohandler.filename = spatialdb/mat_wedge_elastic.spatialdb

observers.observer.data_fields = [displacement, cauchy_stress, cauchy_strain]
observers.observer.trigger.num_skip = 1

auxiliary_subfields.density.basis_order = 0
bulk_rheology.auxiliary_subfields.bulk_modulus.basis_order = 0
bulk_rheology.auxiliary_subfields.shear_modulus.basis_order = 0

derived_subfields.cauchy_strain.basis_order = 1
derived_subfields.cauchy_stress.basis_order = 1


# Mantle
[pylithapp.problem.materials.mantle]
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.label = Properties for mantle
db_auxiliary_field.iohandler.filename = spatialdb/mat_mantle_elastic.spatialdb

observers.observer.data_fields = [displacement, cauchy_stress, cauchy_strain]
observers.observer.trigger.num_skip = 1

auxiliary_subfields.density.basis_order = 0
bulk_rheology.auxiliary_subfields.bulk_modulus.basis_order = 0
bulk_rheology.auxiliary_subfields.shear_modulus.basis_order = 0

derived_subfields.cauchy_strain.basis_order = 1
derived_subfields.cauchy_stress.basis_order = 1


# Crust
[pylithapp.problem.materials.crust]
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.label = Properties for continental crust
db_auxiliary_field.iohandler.filename = spatialdb/mat_crust_elastic.spatialdb

observers.observer.data_fields = [displacement, cauchy_stress, cauchy_strain]
observers.observer.trigger.num_skip = 1

auxiliary_subfields.density.basis_order = 0
bulk_rheology.auxiliary_subfields.bulk_modulus.basis_order = 0
bulk_rheology.auxiliary_subfields.shear_modulus.basis_order = 0

derived_subfields.cauchy_strain.basis_order = 1
derived_subfields.cauchy_stress.basis_order = 1
```

We specify different elastic properties for each material (slab, wedge, mantle, and crust) using `SimpleDB` spatial databases with a single point to specify uniform properties within a material.
We choose `SimpleDB` rather than `UniformDB`, because we will reuse some of these spatial databases for the elastic properties when we use the linear Maxwell viscoelastic constitutive model.

The remaining parameters in the `step01_axialdisp.cfg` file are mostly associated with setting filenames for output of all of the parameters used and version information in a JSON file (`output/step01_axisldisp-parameters.json`) and reporting the progress of the simulation and estimated time of completion (`output/step01_axisldisp-progress.txt`).

```{code-block} console
---
caption: Run Step 1 simulation
---
$ pylith step01_axialdisp.cfg mat_elastic.cfg
```

The simulation will produce ten pairs of HDF5/Xdmf files in the `output` directory:

**step01_axialdisp-domain.h5[.xmf]**  Time series of the solution field over the domain.

**step01_axialdisp-groundsurf.h5[.xmf]**  Time series of the solution field over the ground surface.

**step01_axialdisp-slab_info.h5[.xmf]**  Properties for the slab material.

**step01_axialdisp-slab.h5[.xmf]**  Time series of the state variables (stress and strain) for the slab material.

**step01_axialdisp-wedge)info.h5[.xmf]**  Properties for the wedge material.

**step01_axialdisp-wedge.h5[.xmf]** Time series of the state variables (stress and strain) for the wedge material.

**step01_axialdisp-crust_info.h5[.xmf]**  Properties for the crust material.

**step01_axialdisp-crust.h5[.xmf]** Time series of the tate variables (stress and strain) for the crust material.

**step01_axialdisp-mantle_info.h5[.xmf]** Properties for the mantle material.

**step01_axialdisp-mantle.h5[.xmf]** Time series of the state variables (stress and strain) for the mantle material.

The HDF5 files contain the data and the Xdmf files contain the metadata required by ParaView and Visit (and other visualization tools that use Xdmf files) to access the mesh and data sets in the HDF5 files.

{numref}`fig:example:subduction:3d:step01`, which was created using the ParaView Python script `plot_dispvec.py` (see {ref}`sec:ParaView:Python:scripts` for how to run ParaView Python scripts), displays the magnitude of the displacement field arrows showing the direction and magnitude of the deformation.
Material properties with a positive Poisson's ratio result in vertical deformation along with the axial compression.
The variations in material properties among the properties result in local spatial variations that are most evident in the horizontal displacement components.

:::{figure-md} fig:example:subduction:3d:step01
<img src="figs/subduction3d_step01_soln.*" alt="Solution over the domain for Step 1. The colors indicate the magnitude of the displacement and the arrows indicate the direction with the length of each arrow equal to 10,000 times the magnitude of the displacement." width="100%"/>

Solution over the domain for Step 1. The colors indicate the magnitude of the displacement and the arrows indicate the direction with the length of each arrow equal to 10,000 times the magnitude of the displacement.
:::

## Exercises

* Run the simulation again using multiple cores via the `--nodes=NCORES` argument, replacing `NCORES` with 2 or up to the number of cores on your machine. Examine the PETSc log summary for the various runs to see how the time spent at varies stages changes with the number of cores. Make a plot of runtime versus the number of cores.
* Adjust the material properties in the spatial databases so that the slab is stiffer and the wedge is more compliant. What happens to the solution if you make the materials nearly incompressible? Does this also affect the rate of convergence of the linear solve?
* Change the Dirichlet boundary conditions to impose pure shear instead of axial compression. Hint: You will need to change the boundary conditions on the east, west, north, and south boundaries.
