# Step 2: Prescribed Coseismic Slip and Postseismic Relaxation

In this example we model the postseismic relaxation of the deep slab and mantle resulting from coseismic slip on a fault patch in the central portion of the subduction (top of the slab) interface.
For simplicity we will prescribed uniform slip on the fault patch and use linear Maxwell viscoelastic constitutive models for the slab and mantle.
As the lateral and bottom boundaries are far from the earthquake source, we use roller boundary conditions on these boundaries. We do not expect significant relaxation of stresses on the shallow part of the slab, so we impose a depth-dependent viscosity.
{numref}`fig:example:subduction:3d:step02:diagram` summarizes the problem description.

:::{figure-md} fig:example:subduction:3d:step02:diagram
<img src="figs/subduction3d_step02_diagram.*" alt="Diagram of Step 2 - Prescribed coseismic slip and postseismic relaxation. This quasistatic simulation prescribes uniform slip on the central rupture patch on the subduction interface, depth-dependent viscoelastic relaxation in the slab and mantle, and roller boundary conditions on the lateral (north, south, east, and west) and bottom boundaries." width="100%"/>

Diagram of Step 2 - Prescribed coseismic slip and postseismic relaxation. This quasistatic simulation prescribes uniform slip on the central rupture patch on the subduction interface, depth-dependent viscoelastic relaxation in the slab and mantle, and roller boundary conditions on the lateral (north, south, east, and west) and bottom boundaries.
:::

The `pylithapp.cfg` completely specifies the Dirichlet roller boundary conditions on the five boundaries, so we do not include any boundary condition information in `step02_coseismic.cfg`.
 As discussed in {ref}`sec:example:subduction:3d:organization`, we bundle the parameters for specification of an elastic crust and wedge and viscoelastic slab and mantle in `mat_viscoelastic.cfg`.

We describe the properties of the linear, isotropic Maxwell viscoelastic constitutive model using viscosity in addition to the Vp, Vs, and density used to describe purely linear, isotropic elastic models.
Rather than create a database with all four of these parameters, we leverage the `SimpleDB` spatial databases used by `mat_elastic.cfg` for the elastic properties and simply create a single new spatial database with the depth-dependent viscosity for the slab and mantle.
We use the `CompositeDB` spatial database to combine these two spatial databases into a single spatial database with the material properties.
Rather than using a `SimpleDB` for the depth-dependent viscosity, we use a `SimpleGridDB` spatial database (`spatialdb/mat_viscosity.spatialdb`), which provides faster interpolation using a bilinear search algorithm along each coordinate direction. We use a very large viscosity at depths above 20 km to give behavior that is essentially elastic and decrease it so the Maxwell relaxation time (viscosity divided by the shear modulus) is approximately 200 years at a depth of 30 km, 100 years at a depth of 100 km, and 50 years at a depth of 400 km.
Using linear interpolation results in a piecewise linear variation in the viscosity with depth.

:::{tip}
The `SimpleGridDB` should be used whenever the points in a spatial database can be described with a logically rectangular grid. The grid points along each direction do not need to be uniformly spaced.
:::

In setting the parameters for the `CompositeDB` in `mat_viscoelastic.cfg`, we specify which properties are contained in each of the two spatial databases in the composite database and the type and parameters for each of those spatial databases.
For the slab we have:

```{code-block} cfg
---
caption: Excerpt from `mat_viscoelastic.cfg`
---
[pylithapp.problem.materials.slab]
db_auxiliary_field = spatialdata.spatialdb.CompositeDB
db_auxiliary_field.label = Composite spatial database for slab material properties

observers.observer.data_fields = [displacement, cauchy_stress, cauchy_strain, viscous_strain]
observers.observer.trigger.num_skip = 1

auxiliary_subfields.density.basis_order = 0

derived_subfields.cauchy_stress.basis_order = 1
derived_subfields.cauchy_strain.basis_order = 1

[pylithapp.timedependent.materials.slab.db_auxiliary_field]
# Elastic properties
values_A = [density, vs, vp]
db_A = spatialdata.spatialdb.SimpleDB
db_A.label = Elastic properties for slab
db_A.iohandler.filename = spatialdb/mat_slab_elastic.spatialdb

# Viscoelastic properties
values_B = [viscosity, viscous_strain_xx, viscous_strain_yy, viscous_strain_zz, viscous_strain_xy, viscous_strain_yz, viscous_strain_xz, total_strain_xx, total_strain_yy, total_strain_zz, total_strain_xy, total_strain_yz, total_strain_xz]
db_B = spatialdata.spatialdb.SimpleGridDB
db_B.label = Linear Maxwell viscoelastic properties for slab
db_B.filename = spatialdb/mat_viscosity.spatialdb
db_B.query_type = linear
```

In the simulation specific parameter file `step02_coseismic.cfg`, we specify the parameters for the quasistatic time stepping, the subfields in the solution field, and the coseismic rupture.
By default, PyLith will use implicit time stepping with uniform time steps, so we need only specify the duration and time step size.

```{code-block} cfg
---
caption: Excerpt from `step02_coseismic.cfg` showing the time stepping settings.
---
[pylithapp.problem.formulation.time_step]
initial_dt = 10.0*year
start_time = -10.0*year
end_time = 200.0*year
```

In prescribing coseismic slip on the single fault patch, we create an array with one fault interface and then set its parameters.
Because the edges of the central fault patch are buried within the domain, we need to specify the nodeset that corresponds to the buried edges as well as the nodeset for the entire fault surface.
This ensures that PyLith inserts the cohesive cells and properly terminates the fault surface at the edges.
Just as we do for the boundary conditions and materials, we create an array of components (in this case an array with one fault interface, *slab*, and then refer to those components by name, `pylithapp.problem.interfaces.slab`.
We must also set the discretization information for the fault.

```{code-block} cfg
---
caption: Excerpt from `step02_coseismic.cfg` shwoing the fault settings.
---
[pylithapp.problem]
# We prescribe slip on the slab fault patch.
interfaces = [slab_top]

[pylithapp.problem.interfaces]
slab_top = pylith.faults.FaultCohesiveKin

[pylithapp.problem.interfaces.slab_top]
id = 100
label = fault_slabtop_patch
edge = fault_slabtop_patch_edge

observers.observer.data_fields = [slip]

# The slip time and final slip are defined in spatial databases.
[pylithapp.problem.interfaces.slab_top.eq_ruptures.rupture]
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.label = Coseismic slip
db_auxiliary_field.values = [final_slip_left_lateral, final_slip_reverse, final_slip_opening, initiation_time]
db_auxiliary_field.data = [-1.0*m, 4.0*m, 0.0*m, 9.999*year ]
```

Prescribing the coseismic slip distribution on the fault involves specifying an origin time for the rupture (default is 0.0), and a slip time function along with its parameters.
In this case, we treat the earthquake rupture as just the coseismic slip happening in one time step, so we use a step function for the slip time function (which is the default).
The parameters include the final slip and slip initiation time.
This slip initiation time is relative to the earthquake source origin time, which is 0 by default.
Thus, to specify the time of the slip for a step function, we can either specify the origin time or the slip initiation time; in this case, we use the slip initiation time.
To impose slip beginning at 10 years, we specify a slip initiation time slightly less than 10 years to account for roundoff errors in the nondimensionalization of time.
In Step 4 we will use the origin time.
Because we want uniform slip and a uniform rise time, we use a `UniformDB` spatial database.
Note that we specify oblique slip with 1.0 m of right-lateral motion and 4.0 m of reverse motion.

```{code-block} console
---
caption: Run Step 2 simulation
---
$ pylith step02_coseismic.cfg mat_viscoelastic.cfg solver_fieldsplit.cfg
```
In addition to the ten pairs of HDF5/Xdmf files analogous to those produced in Step 1, we also have two pairs of HDF5/Xdmf files associated with the fault:

**step02_coseismic-domain.h5[.xmf]** Time series of the solution field over the domain.

**step02_coseismic-groundsurf.h5[.xmf]** Time series of the solution field over the ground surface.

**step02_coseismic-slab_info.h5[.xmf]** Properties for the slab material.

**step02_coseismic-slab.h5[.xmf]** Time series of the state variables (stress and strain) for the slab material.

**step02_coseismic-wedge_info.h5[.xmf]** Properties for the wedge material.

**step02_coseismic-wedge.h5[.xmf]** Time series of the state variables (stress and strain) for the wedge material.

**step02_coseismic-crust_info.h5[.xmf]** Properties for the crust material.

**step02_coseismic-crust.h5[.xmf]** Time series of the tate variables (stress and strain) for the crust material.

**step02_coseismic-mantle_info.h5[.xmf]** Properties for the mantle material.

**step02_coseismic-mantle.h5[.xmf]** Time series of the state variables (stress and strain) for the mantle material.

**step02_coseismic-fault-slab_info.h5[.xmf]** Fault orientation and rupture information.

**step02_coseismic-fault-slab.h5[.xmf]** Time series of slip and traction changes.

{numref}`fig:example:subduction:3d:step02`, which was created using the ParaView Python script `plot_dispwarp.py`, displays the magnitude of the displacement field exaggerated by a factor of 10,000 at the final time step (200 yr).
The shallow fault results in deformation that is localized over a small region.

:::{figure-md} fig:example:subduction:3d:step02
<img src="figs/subduction3d_step02_soln.*" alt="Solution over the domain for Step 2 at {math}`t = 200`yr. The colors indicate the magnitude of the displacement and we have exaggerated the deformation by a factor of 10,000." width="100%"/>

Solution over the domain for Step 2 at {math}`t = 200`yr. The colors indicate the magnitude of the displacement and we have exaggerated the deformation by a factor of 10,000.
:::

## Exercises

+ Change the slip from the subduction interface rupture patch to the splay fault rupture patch. Hint: Identify the nodesets for the splay fault patch.
+ Create simultaneous rupture on the subduction interface rupture patch and the splay fault rupture patch.
+ Prescribe coseismic slip on the central patch for splay fault and the subduction interface below the intersection with the splay fault.
  + Implement this without changing any of the nodesets in CUBIT/Trelis. Hint: you will need to create two fault interfaces. What do you notice about the slip at the intersection between the splay fault and slab?
  + Add nodesets in CUBIT/Trelis to create a uniform coseismic slip distribution across the splay fault and on the subduction interface below the splay fault.
