# H5 output schema

This document describes schema version 2 written by
`scripts/simulations/slit_pore/npt.lua`.

## Notation

- `@name` denotes an HDF5 attribute; all other named leaves are datasets.
- `[A, B]` denotes an array shape and `scalar` a scalar dataset or attribute.
- `T_*` is the number of samples at one sampling cadence. The time axes are
  independent: thermodynamics, barostat diagnostics, particle trajectories,
  and structure analysis do not need to contain the same steps.
- `N` is the number of fluid particles, `M` the number of slabs, `Q_o` the
  number of wavevectors for orientation `o`, and `K_o` the corresponding
  number of wavenumber shells. The spatial dimension is `D = 3`.
- A group written as `{step, time, value}` is an H5MD time-dependent element.
  Its `step` and `time` datasets have shape `[T_*]`. Elements produced by the
  same writer may share these datasets through HDF5 hard links.
- Branches marked **optional** are omitted when their sampling interval is
  zero or their required simulation feature is disabled.

## Hierarchy

```text
/
├── h5md/
│   ├── @version                                      [2]
│   ├── author/
│   │   └── @name                                     scalar string
│   ├── creator/
│   │   ├── @name                                     scalar string
│   │   └── @version                                  scalar string
│   └── modules/
│       └── thermodynamics/                           [optional]
│           └── @version                              [2]
│
├── parameters/
│   ├── run/
│   │   ├── @schema_name                              scalar string
│   │   ├── @schema_version                           scalar integer (= 2)
│   │   ├── @experiment_type                          scalar string
│   │   ├── @experiment_name                          scalar string
│   │   ├── @simulation_script                        scalar string
│   │   ├── @output                                   scalar string
│   │   └── @created_at                               scalar string [optional]
│   └── provenance/
│       ├── @repo_commit                              scalar string [optional]
│       ├── @repo_branch                              scalar string [optional]
│       ├── @repo_dirty                               scalar string [optional]
│       ├── @backend                                  scalar string [optional]
│       └── @docker_image                             scalar string [optional]
│
├── geometry/
│   ├── box/
│   │   ├── @dimension                                scalar integer (= 3)
│   │   ├── @boundary                                 [3] strings
│   │   └── edges                                     [3, 3]
│   ├── confinement/
│   │   ├── @axis                                     scalar string
│   │   ├── @pore_width                               scalar
│   │   ├── @analysis_lower                           scalar
│   │   ├── @analysis_upper                           scalar
│   │   ├── @analysis_width                           scalar
│   │   ├── @top_wall_thickness                       scalar
│   │   ├── @bottom_wall_thickness                    scalar
│   │   ├── @flat_wall_padding                        scalar
│   │   ├── @particle_walls_enabled                   scalar integer (0 or 1)
│   │   └── @flat_walls_enabled                       scalar integer (0 or 1)
│   └── slabs/                                        [optional]
│       └── <axis>/
│           ├── @axis                                 scalar string
│           ├── @count                                scalar integer (= M)
│           ├── @requested_width                      scalar
│           ├── @lower                                scalar
│           ├── @upper                                scalar
│           ├── @lower_edge                           [M]
│           ├── @upper_edge                           [M]
│           ├── @center                               [M]
│           └── @volume                               [M]
│
├── particles/                                        [optional]
│   └── fluid/
│       ├── box/                                      automatic H5MD box
│       │   ├── @dimension                            scalar integer (= 3)
│       │   ├── @boundary                             [3] strings
│       │   └── edges                                 [3, 3]
│       ├── position/
│       │   ├── step                                  [T_particles]
│       │   ├── time                                  [T_particles]
│       │   └── value                                 [T_particles, N, 3]
│       └── velocity/
│           ├── step                                  [T_particles]
│           ├── time                                  [T_particles]
│           └── value                                 [T_particles, N, 3]
│
└── observables/
    ├── barostat/                                     [optional]
    │   ├── @target_pressure                          scalar
    │   ├── @lateral_area                             scalar
    │   ├── @maximum_pore_width                       scalar
    │   ├── @piston_mass                              scalar
    │   ├── @damping                                  scalar
    │   ├── @barostat_interval                        scalar integer
    │   ├── @sampling_interval                        scalar integer
    │   ├── @particle_wall_weight                     scalar
    │   ├── @flat_wall_weight                         scalar
    │   ├── @controller_state_fields                  [3] strings
    │   ├── @pressure_component_fields                [3] strings
    │   ├── box/                                      automatic H5MD box
    │   │   ├── @dimension                            scalar integer (= 3)
    │   │   ├── @boundary                             [3] strings
    │   │   └── edges                                 [3, 3]
    │   ├── controller_state/
    │   │   ├── step                                  [T_barostat]
    │   │   ├── time                                  [T_barostat]
    │   │   └── value                                 [T_barostat, 1, 3]
    │   └── pressure_components/
    │       ├── step                                  [T_barostat]
    │       ├── time                                  [T_barostat]
    │       └── value                                 [T_barostat, 1, 3]
    │
    ├── global/
    │   ├── thermodynamics/                           [optional]
    │   │   ├── @dimension                            scalar integer (= 3)
    │   │   ├── particle_number                       scalar
    │   │   ├── density                               scalar
    │   │   ├── volume                                scalar
    │   │   ├── potential_energy/                     {step, time, value[T_thermo]}
    │   │   ├── pressure/                             {step, time, value[T_thermo]}
    │   │   ├── temperature/                          {step, time, value[T_thermo]}
    │   │   ├── virial/                               {step, time, value[T_thermo]}
    │   │   ├── center_of_mass_velocity/              {step, time, value[T_thermo, 3]}
    │   │   └── stress_tensor/                        {step, time, value[T_thermo, 6]}
    │   └── structure/                                [optional]
    │       ├── parallel/
    │       │   ├── density_mode/                     density-mode element
    │       │   └── static_structure_factor/          SSF element
    │       └── normal/
    │           ├── density_mode/                     density-mode element
    │           └── static_structure_factor/          SSF element
    │
    └── slabs/                                        [optional]
        └── <axis>/
            ├── slab_0001/
            │   ├── thermodynamics/
            │   │   ├── @dimension                    scalar integer (= 3)
            │   │   ├── volume                        scalar
            │   │   ├── particle_number/              {step, time, value[T_structure]}
            │   │   ├── density/                      {step, time, value[T_structure]}
            │   │   ├── potential_energy/             {step, time, value[T_structure]}
            │   │   ├── pressure/                     {step, time, value[T_structure]}
            │   │   ├── temperature/                  {step, time, value[T_structure]}
            │   │   ├── virial/                       {step, time, value[T_structure]}
            │   │   └── stress_tensor/                {step, time, value[T_structure, 6]}
            │   └── structure/                        [optional]
            │       ├── parallel/
            │       │   ├── density_mode/             density-mode element
            │       │   └── static_structure_factor/  SSF element
            │       └── normal/
            │           ├── density_mode/             density-mode element
            │           └── static_structure_factor/  SSF element
            ├── slab_0002/
            ├── ...
            └── slab_<MMMM>/
```

`slab_<MMMM>` is zero-padded to four digits, for example `slab_0007`.

## Packed barostat diagnostics

HALMD's phase-space writer safely stores the six scalar barostat diagnostics
as two three-component vectors carried by one synthetic particle. The component
names are explicit attributes on `/observables/barostat`:

```text
@controller_state_fields = [
    pore_width,
    pore_width_velocity,
    normal_pressure
]

@pressure_component_fields = [
    particle_normal_pressure,
    flat_normal_pressure,
    pressure_asymmetry
]
```

Therefore, after removing the singleton particle dimension:

```text
controller_state/value[:, 0, i]
    corresponds to @controller_state_fields[i]

pressure_components/value[:, 0, i]
    corresponds to @pressure_component_fields[i]
```

Inactive wall types retain their pressure channel with value zero. The
`pressure_asymmetry` component is top pressure minus bottom pressure. The
automatic `box/` child is required by the phase-space carrier; canonical
simulation geometry remains under `/geometry`.

## Structure elements

Every `density_mode` group has this layout for its orientation `o`:

```text
density_mode/
├── step                                              [T_structure]
├── time                                              [T_structure]
├── wavevector                                        [Q_o, 3]
└── value                                             [T_structure, Q_o, 2]
```

The last dimension of `value` is `[real, imaginary]`.

Every `static_structure_factor` group has this layout:

```text
static_structure_factor/
├── step                                              [T_structure]
├── time                                              [T_structure]
├── wavenumber                                        [K_o]
└── value                                             [T_structure, K_o, 3]
```

The last dimension of `value` is `[shell_mean, error_of_mean,
wavevector_count]`. These are the SSF values returned by HALMD; they are not
reconstructed from the stored density modes.

The global and slab-resolved structure branches use the same element layout.
`parallel` uses wavevectors parallel to the walls, while `normal` uses
wavevectors along the confinement normal. Slab membership is dynamic, while
HALMD receives the slab occupancy at writer construction as its SSF
normalisation. The time-dependent occupancy is stored in the sibling
`thermodynamics/particle_number` element for analysis-time reweighting or
renormalisation. SSF values for a zero-member slab are undefined and may be
non-finite; analysis should mask those samples.

## Sampling controls and optional branches

All sampling intervals are integers measured in MD steps:

| Argument | Output |
| --- | --- |
| `sampling.thermodynamics` | `/observables/global/thermodynamics` |
| `sampling.barostat` | `/observables/barostat` |
| `sampling.particles` | `/particles/fluid` |
| `sampling.structure` | `/geometry/slabs`, `/observables/slabs`, and structure observables |

The thermodynamics, barostat, and particle writers are connected for the
simulation stages in which their underlying objects are active. Barostat
diagnostics consequently cover only the barostat-controlled stage. Structure
and slab writers are connected for production.

The particle writer covers the complete run. It records step 0 and then every
`sampling.particles` steps through equilibration, barostat, transition, NVT,
and production. HALMD writes periodically extended (unwrapped) positions; the
final state is present only when its step falls on the sampling cadence.

`sampling.structure_observables` further selects the structure payload:

- `all` writes global density modes and SSF, slab thermodynamics, and
  slab-resolved density modes and SSF.
- `thermo` writes slab geometry and slab thermodynamics only. It omits
  `/observables/global/structure` and each slab's `structure/` child.

## Geometry conventions

`pore_width` is the inner wall-to-wall width. The slab analysis bounds cover
the outermost active confinement geometry:

- particle walls extend the corresponding bound by their measured top or
  bottom thickness;
- flat walls extend both bounds by `flat_wall_padding`;
- with both wall models active, each bound is the outermost of the two.

The barostat rejects a pore width that would place either outer particle-wall
surface or either flat-wall plane on or beyond the periodic box boundary. With
box length `L_z`, particle-wall thicknesses `t_top` and `t_bottom`, and flat-wall
padding `p`, the respective strict upper limits are
`L_z - 2 max(t_top, t_bottom)` and `L_z - 2p`. With both wall models enabled,
the stricter of those two limits applies.

The selected slab-axis interval is padded so its analysis width is an integer
multiple of `requested_width`. For `slab.axis = z`, this is the wall-inclusive
confinement interval described above; for `x` or `y`, it is the corresponding
lateral interval. The additional length is divided equally between the lower
and upper bounds. If this symmetrically padded window does not fit inside the
simulation box, analysis setup fails and the simulation box must be increased.
Every slab therefore has exactly `requested_width` along the slab axis, and all
slabs have equal volume for a given orthogonal cross-sectional area. Flat walls
have no particle records and remain represented by `/geometry/confinement` and
the barostat pore-width history.
