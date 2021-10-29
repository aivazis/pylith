# Step 8: Stress Field Due to Gravitational Body Forces

This example demonstrates the use of gravitational body forces as well as the use of initial stresses to balance the body forces.
This involves enabling gravity within our domain with Dirichlet roller boundary conditions on the lateral and bottom boundaries; we do not include faults in this example.
We also demonstrate what happens when the initial stresses are not in balance with the gravitational stresses, and show how viscoelastic problems with gravitational stresses will, in general, not reach a steady-state solution.
The example is divided into three sub-problems:

Step 8a
: Gravitational body forces with 3D density variations in elastic materials and initial stresses for a uniform density.

Step 8b
: Gravitational body forces with 3D density variations in incompressible elastic materials.

Step 8c
: Gravitational body forces with 3D density variations in elastic and viscoelastic materials and initial stresses from Step 8b plus finite strain formulation (does not reach a steady-state solution).

## Step 8a

For Step 8a we apply gravitational stresses and attempt to balance these with analytically computed stresses consistent with the density of the mantle.
Since the actual density is not uniform, the stresses are out of balance and we end up with some deformation.
In `step08a_gravity_refstate.cfg` we turn on gravity and set the total time to zero (there is no time dependence in this model).
We set the basis order of the displacement subfield to 1, because we expect little deformation.

:::{important}
Quasistatic problems do not have a well-defined density (inertia) scale.
For this static simulation we adjust the time scale and time step to give a density scale close to unity.
:::



```{code-block} cfg
---
caption: Excerpt from `step08a_gravity_refstate.cfg` showing the gravity settings.
---
[pylithapp.problem]
initial_dt = 1.0*s
start_time = -1.0*s
end_time = 0.0*s
normalizer.relaxation_time = 1.0*s

gravity_field = spatialdata.spatialdb.GravityField
gravity_field.gravity_dir = [0.0, 0.0, -1.0]

defaults.quadrature_order = 1

[pylithapp.problem.solution.subfields.displacement]
basis_order = 1
```

Our initial stress field corresponds to {math}`\sigma_{xx} = \sigma_{yy} = \sigma_{zz} = \rho_{mantle} g z` for all four materials, where {math}`\rho_{mantle}` is the density of the mantle, {math}`g` is the acceleration due to gravity, and {math}`z` is elevation.
With only two control points necessary to describe this linear variation, we use `SimpleDB` spatial databases for all four materials.

```{code-block} cfg
---
caption: Excerpt from `step08a_gravity_refstate.cfg` showing the settings to specify the reference stress.
---
[pylithapp.problem.materials.slab]
db_auxiliary_field.iohandler.filename = spatialdb/mat_slab_elastic_refstate.spatialdb
db_auxiliary_field.query_type = linear

[pylithapp.problem.materials.slab.bulk_rheology]
use_reference_state = True
auxiliary_subfields.reference_stress.basis_order = 1
auxiliary_subfields.reference_strain.basis_order = 0


# Wedge
[pylithapp.problem.materials.wedge]
db_auxiliary_field.iohandler.filename = spatialdb/mat_wedge_elastic_refstate.spatialdb
db_auxiliary_field.query_type = linear

[pylithapp.problem.materials.wedge.bulk_rheology]
use_reference_state = True
auxiliary_subfields.reference_stress.basis_order = 1
auxiliary_subfields.reference_strain.basis_order = 0


# Mantle
[pylithapp.problem.materials.mantle]
db_auxiliary_field.iohandler.filename = spatialdb/mat_mantle_elastic_refstate.spatialdb
db_auxiliary_field.query_type = linear

[pylithapp.problem.materials.mantle.bulk_rheology]
use_reference_state = True
auxiliary_subfields.reference_stress.basis_order = 1
auxiliary_subfields.reference_strain.basis_order = 0


# Crust
[pylithapp.problem.materials.crust]
db_auxiliary_field.iohandler.filename = spatialdb/mat_crust_elastic_refstate.spatialdb
db_auxiliary_field.query_type = linear

[pylithapp.problem.materials.crust.bulk_rheology]
use_reference_state = True
auxiliary_subfields.reference_stress.basis_order = 1
auxiliary_subfields.reference_strain.basis_order = 0
```

```{code-block} console
---
caption: Run Step 8a simulation
---
$ pylith step08a_gravity_refstate.cfg mat_elastic.cfg
```

The simulation will generate HDF5/Xdmf files beginning with `step08a_gravity_refstate`:

**step08a_gravity_refstate-domain.h5[.xmf]**  Time series of the solution field over the domain.

**step08a_gravity_refstate-groundsurf.h5[.xmf]**  Time series of the solution field over the ground surface.

**step08a_gravity_refstate-slab_info.h5[.xmf]**  Properties for the slab material.

**step08a_gravity_refstate-slab.h5[.xmf]**  Time series of the state variables (stress and strain) for the slab material.

**step08a_gravity_refstate-wedge_info.h5[.xmf]**  Properties for the wedge material.

**step08a_gravity_refstate-wedge.h5[.xmf]**  Time series of the state variables (stress and strain) for the wedge material.

**step08a_gravity_refstate-crust_info.h5[.xmf]**  Properties for the crust material.

**step08a_gravity_refstate-crust.h5[.xmf]**  Time series of the tate variables (stress and strain) for the crust material.

**step08a_gravity_refstate-mantle_info.h5[.xmf]**  Properties for the mantle material.

**step08a_gravity_refstate-mantle.h5[.xmf]**  Time series of the state variables (stress and strain) for the mantle material.

The solution shows deformation consistent with the difference between the actual and assumed uniform density for the reference stress.
The slab subsides while the crust undergoes uplift due to the differences in density relative to the mantle.
{numref}`fig:example:subduction:3d:step08a` shows the deformed mesh visualized with the `plot_dispwarp.py` ParaView Python script.

:::{figure-md} fig:example:subduction:3d:step08a
<img src="figs/subduction3d_step08a_soln.*" alt="Solution for Step 8a. The deformation has been exaggerated by a factor of 500 and the colors highlight the vertical displacement component. The crustal material in the east is less dense than the assumed mantle material for initial stresses, while the slab material in the west is more dense. The result is uplift in the east and subsidence in the west." width="100%"/>

Solution for Step 8a.
The deformation has been exaggerated by a factor of 500 and the colors highlight the vertical displacement component.
The crustal material in the east is less dense than the assumed mantle material for initial stresses, while the slab material in the west is more dense.
The result is uplift in the east and subsidence in the west.
:::

## Step 8b

Step 8b is similar to Step 8a, but we replace the elastic material with the incompressible elastic material.
With the incompressible elastic material we do not need to set the reference stress to avoid volumetric deformation.
In general, the incompressible elastic meterial yields a better approximation of the stress state compared to an elastic material with a reference stress because it retains the original geometry better.

We consolidate the settings for the incompressible elastic materials in `mat_elastic_incompressible.cfg`.
We change the materials to `IncompressibleElasticity` and the bulk rheologies to `IsotropicLinearIncompElasticity`.

```{code-block} cfg
---
caption: Excerpt from `mat_elastic_incompressible.cfg` showing the settings for the materials and bulk rheologies.
---
[pylithapp.problem.materials]
slab = pylith.materials.IncompressibleElasticity
wedge = pylith.materials.IncompressibleElasticity
crust = pylith.materials.IncompressibleElasticity
mantle = pylith.materials.IncompressibleElasticity

slab.bulk_rheology = pylith.materials.IsotropicLinearIncompElasticity
wedge.bulk_rheology = pylith.materials.IsotropicLinearIncompElasticity
crust.bulk_rheology = pylith.materials.IsotropicLinearIncompElasticity
mantle.bulk_rheology = pylith.materials.IsotropicLinearIncompElasticity
```

:::{note}
Currently, the only bulk rheology implemented for the incompressible elastic material is isotropic linear elasticity.
:::

We set the properties for the incompressible elastic materials in the same way as we do for the elastic materials.

```{code-block} cfg
---
caption: Excerpt from `mat_elastic_incompressible.cfg` showing the settings for the slap material.
---
[pylithapp.problem.materials.slab]
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.label = Properties for subducting slab
db_auxiliary_field.values = [density, vs, vp]
db_auxiliary_field.data = [3400*kg/m**3, 4.5*km/s, 1.0e+12*km/s]

observers.observer.data_fields = [displacement, cauchy_stress, cauchy_strain]
observers.observer.trigger.num_skip = 1

auxiliary_subfields.density.basis_order = 0
bulk_rheology.auxiliary_subfields.bulk_modulus.basis_order = 0
bulk_rheology.auxiliary_subfields.shear_modulus.basis_order = 0

derived_subfields.cauchy_strain.basis_order = 1
derived_subfields.cauchy_stress.basis_order = 1
```


The formulation for the incompressible elastic material uses displacement and pressure solution subfields.
We expect the pressure to increase nearly linearly with depth so we use a basis order of 1 for the pressure field.
We expect little deformation, so we use a basis order of 1 for the displacement field.

```{code-block} cfg
---
caption: Excerpt from `step08b_gravity_incompressible.cfg` showing the settings for the solution field.
---
[pylithapp.problem]
# We use the predefined container with displacement and pressure (mean
# compressive stress) subfields for the solution field.
solution = pylith.problems.SolnDispPres

defaults.quadrature_order = 1

[pylithapp.problem.solution.subfields]
# We reduce the basis order to 1 because we expect little or no
# deformation with incompressible elasticity.
displacement.basis_order = 1
pressure.basis_order = 1
```

With the addition of the pressure subfield, we add a boundary condition to prescribe zero pressure at the ground surface.

```{code-block} cfg
---
caption: Excerpt from `step08b_gravity_incompressible.cfg` showing the settings for the pressure Dirichlet boundary condition on the ground surface.
---
[pylithapp.problem]
bc = [bc_xneg, bc_xpos, bc_yneg, bc_ypos, bc_zneg, bc_zpos]
bc.bc_zpos = pylith.bc.DirichletTimeDependent

# This BC must be fully specified since it is not included in pylithapp.cfg.
[pylithapp.problem.bc.bc_zpos]
constrained_dof = [0]
label = boundary_zpos
field = pressure
db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.label = Dirichlet BC for pressure on ground surface

auxiliary_subfields.initial_amplitude.basis_order = 0

observers.observer.data_fields = [pressure]
```

We also update the PETSc solver settings to use algebraic multigrid preconditioning on the displacement and pressure separately.

```{code-block} cfg
---
caption: Excerpt from `step08b_gravity_incompressible.cfg` showing the PETSc solver settings.
---
[pylithapp.petsc]
pc_type = fieldsplit
pc_fieldsplit_type = schur
pc_fieldsplit_schur_fact_type = full
pc_fieldsplit_schur_precondition = full
fieldsplit_displacement_pc_type = ml
fieldsplit_pressure_pc_type = ml
```

```{code-block} console
---
caption: Run Step 8b simulation
---
$ pylith step08b_gravity_incompressible.cfg mat_elastic_incompressible.cfg
```

## Step 8c

:::{error}
This example has not been updated for v3 because the small-strain formulation has not been re-implemented in the multiphysics formulation.
::: 


%Because the initial stresses are consistent with the variations in density across the materials, the initial stresses will satisfy equilibrium and there will be essentially no deformation.
%We use the Python script `generate_initial_stress.py`, located in the `spatialdb` directory, to postprocess the output from Step 8a and generate the initial stress spatial database.
%Note that this script uses the Python interface to the spatialdata package to write the spatial database; this is much easier than writing a script to format the data to conform to the format of the spatial database.
%The spatial database will contain the stresses at each vertex of our unstructured mesh, so the points are not on a logical grid and we must use a `SimpleDB`.

%```{code-block} console
%---
%caption: Generate the initial stresses for Step 8b
%---
%# From the examples/3d/subduction directory, change to the spatialdb subdirectory.
%$ cd spatialdb
%$ ./generate_initial_stress.py
%```

%This will create spatial databases containing initial stresses for each of the four materials.

%In the `step08b.cfg` file we specify the *SimpleDB* spatial database for each material (they are now material specific).
%With points at each cell centroid, we use nearest interpolation (default) rather than linear interpolation; this is a small approximation but it is much faster than using linear interpolation in this unstructured set of points.

%```{code-block} cfg
%---
%caption: Excerpt from `step08b.cfg`
%---
%```

%```{code-block} console
%---
%caption: Run Step 8b simulation
%---
%$ pylith step08b.cfg mat_elastic.cfg solver_algebraicmultigrid.cfg
%```
%This simulation will produce files in the `output` directory analogous to Step 8a.

%When we compare the resulting elastic displacements with those of Step 8a, we find that there is essentially no displacement, as seen in {numref}`fig:example:subduction:3d:step08b`.

%:::{figure-md} fig:example:subduction:3d:step08b
%<img src="figs/subduction3d_step08b_soln.*" alt="Solution for Step 8b. In this case the initial stresses satisfy the governing equation, so there is no deformation." width="100%"/>

%Solution for Step 8b. In this case the initial stresses satisfy the governing equation, so there is no deformation.
%:::

%In this example we use linear Maxwell viscoelastic models in place of the elastic models for the slab and mantle.
%We also use the small strain formulation (*ImplicitLgDeform*) so that the deformed configuration is taken into account; Steps 8a and 8b use the default *Implicit* infinitesimal strain formulation.
%The small strain formulation should generally be used for viscoelastic problems with gravity where you need accurate estimates of the vertical deformation.

%:::{warning}
%The shear stress variations in the initial stresses will cause the viscoelastic materials to drive viscous flow, resulting in time-dependent deformation.
%As long as the elastic materials impose deviatoric stresses in the viscoelastic materials through continuity of strain, the viscoelastic materials will continue to flow.
%**As a result, in this case and many other simulations with viscoelastic materials and gravitational body forces, it is difficult to find a steady state solution.**
%:::

%The only difference between the parameters in `step08b.cfg` and `step08c.cfg` is in the formulation setting and the simulation time:

## Exercises

+ What happens in Step 8a if we use a different reference density to compute our initial stresses?
+ For Step 8c what happens if you:
  + Run the simulation for a longer period of time?
  + Change the viscoelastic properties? For example, reduce the viscosity, make all materials viscoelastic, switch to a power-law rheology, etc.
  + Is it possible to find a better initial stress state for Step 8c?
