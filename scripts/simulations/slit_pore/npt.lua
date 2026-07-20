-- Lennard-Jones fluid in a slit pore with optional NPAT pressure control.

local mdsim = halmd.mdsim
local observables = halmd.observables
local readers = halmd.io.readers
local writers = halmd.io.writers
local utility = halmd.utility

package.path = utility.abspath("../../lib") .. "/?.lua;" .. package.path

local lennard_jones = require("lennard_jones")
local slab_analysis = require("slab_analysis")
local ssf_analysis = require("ssf_analysis")
local wavevectors_mod = require("wavevectors")
local h5_metadata = require("h5_metadata")
local data_dir = utility.abspath("../../../data")
local barostat = require("barostat")

local function is_finite(value)
    return type(value) == "number"
        and value == value
        and value < math.huge
        and value > -math.huge
end

local function require_finite(name, value)
    if not is_finite(value) then
        error(("%s must be finite: %s"):format(name, tostring(value)), 3)
    end
end

local function require_positive(name, value)
    require_finite(name, value)
    if value <= 0 then
        error(name .. " must be positive", 3)
    end
end

local function require_nonnegative(name, value)
    require_finite(name, value)
    if value < 0 then
        error(name .. " must not be negative", 3)
    end
end

local function validate_inputs(args)
    for _, entry in ipairs({
        {"pore length", args.pore_length},
        {"initial pore width", args.pore_width},
        {"simulation box width", args.box_width},
        {"fluid density", args.fluid_density},
        {"temperature", args.temperature},
        {"fluid and particle-wall cutoff", args.cutoff},
        {"flat-wall cutoff", args.flat_cutoff},
        {"cutoff smoothing", args.smoothing},
        {"timestep", args.timestep},
        {"slab width", args.slab.width},
        {"piston mass", args.barostat.mass},
        {"maximum pore-width change", args.barostat.max_width_change}
    }) do
        require_positive(entry[1], entry[2])
    end

    for _, entry in ipairs({
        {"flat-wall offset", args.flat_offset},
        {"heat-bath collision rate", args.rate},
        {"initial equilibration time", args.equilibrium_time},
        {"barostat time", args.barostat_time},
        {"NVT equilibration time", args.nvt_equilibrium_time},
        {"production time", args.production_time},
        {"barostat damping", args.barostat.damping},
        {"barostat average time", args.barostat.average_time},
        {"maximum transition time", args.barostat.max_transition_time},
        {"NVT switch tolerance", args.barostat.switch_tolerance},
        {"particle-wall weight", args.weights.particle},
        {"flat-wall weight", args.weights.flat}
    }) do
        require_nonnegative(entry[1], entry[2])
    end

    require_finite("target pressure", args.barostat.target_pressure)
    require_positive("barostat interval", args.barostat.interval)
    if args.barostat.interval % 1 ~= 0 then
        error("barostat interval must be an integer", 2)
    end

    if args.weights.particle == 0 and args.weights.flat == 0 then
        error("particle walls and/or flat walls are required", 2)
    end
end

local function validate_pore_width(name, width, max_width)
    require_positive(name, width)
    if width >= max_width then
        error(("%s (%g) reaches or exceeds maximum pore width (%g)")
            :format(name, width, max_width), 3)
    end
end

local function read_sample(path, location, fields)
    local file = readers.h5md({path = path})
    local reader, sample = observables.phase_space.reader({
        file = file,
        location = location,
        fields = fields
    })
    reader:read_at_step(-1)
    local edges = mdsim.box.reader({
        file = file,
        location = location
    })
    file:close()
    return edges, sample
end

local function clean_wall_particle(particle, label)
    local clean = mdsim.particle({
        dimension = particle.dimension,
        particles = particle.nparticle,
        species = particle.nspecies,
        memory = particle.memory,
        precision = particle.precision,
        label = label
    })

    clean.data.position = particle.data.position
    clean.data.species = particle.data.species

    return clean
end

local function split_wall_particles(box, particle_obstacles)
    local corner = box:lowest_corner()
    local length = box.length
    local center = 0

    local top_bbox = mdsim.geometries.cuboid({
        lowest_corner = {corner[1], corner[2], center},
        length = {
            length[1],
            length[2],
            corner[3] + length[3] - center
        }
    })

    local top_selection = mdsim.particle_groups.region({
        particle = particle_obstacles,
        geometry = top_bbox,
        selection = "included",
        label = "top_wall_selection",
        fluctuating = false
    })

    local bottom_selection = mdsim.particle_groups.region({
        particle = particle_obstacles,
        geometry = top_bbox,
        selection = "excluded",
        label = "bottom_wall_selection",
        fluctuating = false
    })

    assert(top_selection.size > 0, "top wall has no particles")
    assert(bottom_selection.size > 0, "bottom wall has no particles")

    local particle_wall_top = top_selection:to_particle({
        label = "wall_top"
    })
    local particle_wall_bottom = bottom_selection:to_particle({
        label = "wall_bottom"
    })

    -- Work around corrupt double-single low words from GPU to_particle().
    particle_wall_top = clean_wall_particle(
        particle_wall_top,
        "wall_top"
    )
    particle_wall_bottom = clean_wall_particle(
        particle_wall_bottom,
        "wall_bottom"
    )

    local group_wall_top = mdsim.particle_groups.all({
        particle = particle_wall_top,
        label = "wall_top"
    })
    local group_wall_bottom = mdsim.particle_groups.all({
        particle = particle_wall_bottom,
        label = "wall_bottom"
    })

    return particle_wall_top, particle_wall_bottom,
       group_wall_top, group_wall_bottom
end

local function z_bounds(particle)
    local z_min = math.huge
    local z_max = -math.huge

    for _, position in ipairs(particle.data.position) do
        local z = position[3]
        z_min = math.min(z_min, z)
        z_max = math.max(z_max, z)
    end

    return z_min, z_max
end

local function load_wall_particles(args, box)
    local obstacle_edges, obstacle_sample = read_sample(
        args.particle_wall_input,
        {"particles", "obstacles"},
        {"position", "species"}
    )

    local obstacle_box = mdsim.box({
        edges = obstacle_edges
    })

    local obstacle_length = obstacle_box.length
    local simulation_box_length = box.length

    assert(
        args.pore_length < obstacle_length[1]
            and args.pore_length < obstacle_length[2],
        ("pore_length (%g) must be smaller than obstacle box (%g x %g)")
            :format(
                args.pore_length,
                obstacle_length[1],
                obstacle_length[2]
            )
    )
    assert(
        simulation_box_length[3] > obstacle_length[3],
        ("simulation box width (%g) must be larger than obstacle box width (%g)")
            :format(simulation_box_length[3], obstacle_length[3])
    )

    local particle_obstacles_all = mdsim.particle({
        dimension = obstacle_sample.dimension,
        particles = obstacle_sample.nparticle,
        species = obstacle_sample.nspecies,
        label = "obstacles_all"
    })

    local group_obstacles_all = mdsim.particle_groups.all({
        particle = particle_obstacles_all,
        label = "obstacles_all"
    })

    observables.phase_space({
        box = obstacle_box,
        group = group_obstacles_all
    }):set(obstacle_sample)

    local lateral_bbox = mdsim.geometries.cuboid({
        lowest_corner = {
            -0.5 * simulation_box_length[1],
            -0.5 * simulation_box_length[2],
            -0.5 * simulation_box_length[3]
        },
        length = {
            simulation_box_length[1],
            simulation_box_length[2],
            simulation_box_length[3]
        }
    })

    local lateral_selection = mdsim.particle_groups.region({
        particle = particle_obstacles_all,
        geometry = lateral_bbox,
        selection = "included",
        label = "obstacles_lateral_selection",
        fluctuating = false
    })

    assert(
        lateral_selection.size > 0,
        "no obstacle particles are chosen after lateral crop"
    )

    local particle_obstacles = lateral_selection:to_particle({
        label = "obstacles"
    })

    local particle_wall_top,
    particle_wall_bottom,
    group_wall_top,
    group_wall_bottom =
        split_wall_particles(box, particle_obstacles)

    local bottom_outer, bottom_inner = z_bounds(particle_wall_bottom)
    local top_inner, top_outer = z_bounds(particle_wall_top)

    local top_thickness = top_outer - top_inner
    local bottom_thickness = bottom_inner - bottom_outer

    return {
        top = {
            particle = particle_wall_top,
            group = group_wall_top,
            thickness = top_thickness,
            inner = top_inner
        },
        bottom = {
            particle = particle_wall_bottom,
            group = group_wall_bottom,
            thickness = bottom_thickness,
            inner = bottom_inner
        }
    }
end

local function position_particle_walls(particle_walls, box, width)
    local half_width = 0.5 * width
    local top = particle_walls.top
    local bottom = particle_walls.bottom

    top.particle:shift_position(
        box,
        {0, 0, half_width - top.inner}
    )
    bottom.particle:shift_position(
        box,
        {0, 0, -half_width - bottom.inner}
    )

    top.inner = half_width
    bottom.inner = -half_width
end

local function maximum_pore_width(
    args,
    box,
    particle_walls,
    flat_walls_enabled
)
    local max_width = math.huge

    if particle_walls then
        max_width = math.min(
            max_width,
            box.length[3] - 2 * math.max(
                particle_walls.top.thickness,
                particle_walls.bottom.thickness
            )
        )
    end

    if flat_walls_enabled then
        max_width = math.min(
            max_width,
            box.length[3] - 2 * args.flat_offset
        )
    end

    require_positive("maximum pore width", max_width)
    return max_width
end

local function analysis_region(args, width, particle_walls, flat_walls, box)
    local half_width = 0.5 * width
    local lower = -half_width
    local upper = half_width
    local top_wall_thickness = 0
    local bottom_wall_thickness = 0
    local flat_wall_padding = 0

    if particle_walls then
        top_wall_thickness = particle_walls.top.thickness
        bottom_wall_thickness = particle_walls.bottom.thickness
        lower = lower - bottom_wall_thickness
        upper = upper + top_wall_thickness
    end

    if flat_walls then
        flat_wall_padding = args.flat_offset
        lower = math.min(lower, -half_width - flat_wall_padding)
        upper = math.max(upper, half_width + flat_wall_padding)
    end

    local box_corner = box:lowest_corner()
    local box_lower = box_corner[3]
    local box_upper = box_lower + box.length[3]
    assert(
        lower >= box_lower and upper <= box_upper,
        ("analysis range [%g, %g] exceeds simulation box [%g, %g]")
            :format(lower, upper, box_lower, box_upper)
    )

    return {
        lowest_corner = {
            -0.5 * args.pore_length,
            -0.5 * args.pore_length,
            lower
        },
        length = {
            args.pore_length,
            args.pore_length,
            upper - lower
        }
    }, {
        axis = "z",
        pore_width = width,
        analysis_lower = lower,
        analysis_upper = upper,
        top_wall_thickness = top_wall_thickness,
        bottom_wall_thickness = bottom_wall_thickness,
        flat_wall_padding = flat_wall_padding,
        particle_walls_enabled = particle_walls ~= nil,
        flat_walls_enabled = flat_walls ~= nil
    }
end

local function create_particle_wall_potential(args, weight)
    local particle_wall_potential = mdsim.potentials.pair.lennard_jones({
        epsilon = {
            {0, weight},
            {weight, 0}
        },
        sigma = {
            {1, 1},
            {1, 1}
        },
        label = "particle wall potential"
    })

    particle_wall_potential = particle_wall_potential:truncate({
        "smooth_r4",
        cutoff = args.cutoff,
        h = args.smoothing
    })

    return particle_wall_potential
end

local function create_flat_wall_potential(args, side, offset, weight)
    local normals = {
        top = {0, 0, 1},
        bottom = {0, 0, -1}
    }

    return mdsim.potentials.external.planar_wall({
        surface_normal = normals[side],
        offset = {offset},
        sigma = 1,
        epsilon = weight,
        wetting = 0,
        cutoff = args.flat_cutoff,
        smoothing = args.smoothing,
        species = 1
    })
end

function main(args)
    validate_inputs(args)

    local weights = args.weights
    local particle_walls_enabled = weights.particle > 0
    local flat_walls_enabled = weights.flat > 0
    local barostat_enabled = args.barostat_time > 0
    local structure_observables = args.sampling.structure_observables

    if structure_observables ~= "all"
        and structure_observables ~= "thermo"
    then
        error("sampling.structure_observables must be 'all' or 'thermo'")
    end

    for _, name in ipairs({
        "thermodynamics",
        "barostat",
        "particles",
        "structure"
    }) do
        if args.sampling[name] < 0 then
            error("sampling." .. name .. " must not be negative")
        end
    end

    local piston
    local lateral_area = args.pore_length * args.pore_length
    local pore_volume = function()
        local width = piston and piston.width or args.pore_width
        return lateral_area * width
    end

    local nparticle = math.ceil(pore_volume() * args.fluid_density)

    local box = mdsim.box({
        length = {
            args.pore_length,
            args.pore_length,
            args.box_width
        }
    })

    local particle_fluid = mdsim.particle({
        particles = nparticle,
        dimension = 3,
        species = 1,
        label = "fluid"
    })
    local group_fluid = mdsim.particle_groups.all({
        particle = particle_fluid,
        label = "fluid"
    })

    -- Create and initialize fluid particles.
    mdsim.positions.lattice({
        box = box,
        particle = particle_fluid,
        slab = {
            0.95,
            0.95,
            0.90 * args.pore_width / args.box_width
        }
    }):set()
    mdsim.velocities.boltzmann({
        particle = particle_fluid,
        temperature = args.temperature
    }):set()
    particle_fluid.data["species"] = utility.repeat_element(
        0,
        particle_fluid.nparticle
    )

    local particle_walls
    if particle_walls_enabled then
        readers.h5md.check(args.particle_wall_input)
        particle_walls = load_wall_particles(args, box)
    end

    local max_width = maximum_pore_width(
        args,
        box,
        particle_walls,
        flat_walls_enabled
    )
    validate_pore_width("initial pore width", args.pore_width, max_width)

    if particle_walls then
        position_particle_walls(particle_walls, box, args.pore_width)
    end

    -- Register wall forces first. Their order is recorded for recovering
    -- individual wall contributions from cumulative fluid-force snapshots.
    local force_order = {}

    -- Create particle-wall forces.
    local particle_wall_potential
    if particle_walls_enabled then
        particle_wall_potential = create_particle_wall_potential(
            args,
            weights.particle
        )

        local force_top_particle = mdsim.forces.pair({
            box = box,
            particle = {particle_fluid, particle_walls.top.particle},
            potential = particle_wall_potential
        })

        force_order[#force_order + 1] = {
            label = "top_particle",
            force = force_top_particle
        }

        local force_bottom_particle = mdsim.forces.pair({
            box = box,
            particle = {particle_fluid, particle_walls.bottom.particle},
            potential = particle_wall_potential
        })

        force_order[#force_order + 1] = {
            label = "bottom_particle",
            force = force_bottom_particle
        }
    end

    -- Create flat-wall forces.
    local flat_walls
    if flat_walls_enabled then
        local flat_wall_offset = 0.5 * args.pore_width + args.flat_offset

        local top_flat_potential = create_flat_wall_potential(
            args, "top", flat_wall_offset, weights.flat
        )

        local top_flat_force = mdsim.forces.external({
            box = box,
            particle = particle_fluid,
            potential = top_flat_potential
        })

        force_order[#force_order + 1] = {
            label = "top_flat",
            force = top_flat_force
        }

        local bottom_flat_potential = create_flat_wall_potential(
            args, "bottom", flat_wall_offset, weights.flat
        )

        local bottom_flat_force = mdsim.forces.external({
            box = box,
            particle = particle_fluid,
            potential = bottom_flat_potential
        })

        force_order[#force_order + 1] = {
            label = "bottom_flat",
            force = bottom_flat_force
        }

        flat_walls = {
            top = top_flat_potential,
            bottom = bottom_flat_potential
        }
    end

    -- Register fluid-fluid last. It is part of the physical force evaluation
    -- but intentionally excluded from wall-force recovery.
    local force_fluid_fluid = lennard_jones.create_pair_force({
        box = box,
        particle = particle_fluid,
        cutoff = args.cutoff,
        smoothing = args.smoothing
    })

    local fluid_thermodynamics = observables.thermodynamics({
        box = box,
        group = group_fluid,
        volume = pore_volume
    })

    -- Add the initial velocity-Verlet NVT integrator.
    local equilibration_integrator = mdsim.integrators.verlet_nvt_boltzmann({
        box = box,
        particle = particle_fluid,
        timestep = args.timestep,
        temperature = args.temperature,
        rate = args.rate * 10
    })

    -- H5MD file writer
    local file = writers.h5md({
        path = ("%s.h5"):format(args.output),
        overwrite = args.overwrite
    })
    h5_metadata.write_run(file, args, {
        schema_version = 2,
        experiment_type = "npat",
        simulation_script = "scripts/simulations/slit_pore/npt.lua"
    })
    box:writer({file = file, location = {"geometry"}})

    local equilibrium_steps = math.ceil(args.equilibrium_time / args.timestep)
    local barostat_steps = math.ceil(args.barostat_time / args.timestep)
    local nvt_equilibrium_steps = math.ceil(args.nvt_equilibrium_time / args.timestep)
    local production_steps = math.ceil(args.production_time / args.timestep)
    local transition_to_nvt = args.barostat.transition_to_nvt
    assert(
        not transition_to_nvt or barostat_enabled,
        "barostat-to-NVT transition requires positive barostat time"
    )
    local max_transition_steps = transition_to_nvt
        and math.ceil(
            args.barostat.max_transition_time / args.timestep
        )
        or 0

    local total_steps =
        equilibrium_steps
        + barostat_steps
        + max_transition_steps
        + nvt_equilibrium_steps
        + production_steps

    if args.sampling.thermodynamics > 0 then
        fluid_thermodynamics:writer({
            file = file,
            location = {"observables", "global", "thermodynamics"},
            fields = {
                "potential_energy",
                "pressure",
                "temperature",
                "center_of_mass_velocity",
                "stress_tensor",
                "virial"
            },
            every = args.sampling.thermodynamics
        })
    end

    if args.sampling.particles > 0 then
        observables.phase_space({
            box = box,
            group = group_fluid
        }):writer({
            file = file,
            location = {"particles", "fluid"},
            fields = {"position", "velocity"},
            every = args.sampling.particles
        })
    end

    -- Sample the initial state.
    observables.sampler:sample()
    -- Estimate the remaining runtime.
    observables.runtime_estimate({steps = total_steps})

    -- Run the initial fixed-width NVT equilibration.
    observables.sampler:run(equilibrium_steps)
    equilibration_integrator:disconnect()

    -- Add the velocity-Verlet NVT integrator used by later stages.
    local barostat_integrator = mdsim.integrators.verlet_nvt_boltzmann({
        box = box,
        particle = particle_fluid,
        timestep = args.timestep,
        temperature = args.temperature,
        rate = args.rate
    })

    local average_window_steps = 0
    if transition_to_nvt then
        average_window_steps = math.ceil(
            args.barostat.average_time / args.timestep
        )
        assert(
            average_window_steps <= barostat_steps,
            "barostat average time must not exceed barostat time"
        )
        assert(
            average_window_steps >= args.barostat.interval,
            "barostat average time must span at least one update"
        )
    end

    local width_sum = 0
    local width_samples = 0

    local function collect_average_width(pressure_barostat, step)
        if step > barostat_steps - average_window_steps
            and step <= barostat_steps
        then
            width_sum = width_sum + pressure_barostat.width
            width_samples = width_samples + 1
        end
    end

    if barostat_enabled then
        piston = barostat.create({
            box = box,
            particle = particle_fluid,
            thermodynamics = fluid_thermodynamics,
            wall_forces = force_order,
            particle_walls = particle_walls,
            flat_walls = flat_walls,
            flat_offset = args.flat_offset,
            area = lateral_area,
            initial_width = args.pore_width,
            max_width = max_width,
            mass = args.barostat.mass,
            target_pressure = args.barostat.target_pressure,
            damping = args.barostat.damping,
            max_width_change = args.barostat.max_width_change,
            timestep = args.timestep,
            interval = args.barostat.interval,
            on_update = transition_to_nvt
                and collect_average_width
                or nil,
            diagnostics = {
                file = file,
                interval = args.sampling.barostat,
                particle_wall_weight = weights.particle,
                flat_wall_weight = weights.flat
            }
        })
        piston:connect()
    end

    -- Run the NPAT equilibration stage.
    observables.sampler:run(barostat_steps)

    -- Run the optional NPAT-to-NVT transition.
    if transition_to_nvt then
        assert(width_samples > 0, "no pore-width samples were collected")
        local average_width = width_sum / width_samples
        local transition_steps = 0

        while math.abs(piston.width - average_width)
                > args.barostat.switch_tolerance
            do
                if transition_steps + args.barostat.interval > max_transition_steps then
                    error("barostat-to-NVT transition tolerance was not reached")
                end

                observables.sampler:run(args.barostat.interval)
                transition_steps = transition_steps + args.barostat.interval
            end
        end

    if piston then
        piston:disconnect()
    end

    -- Run the fixed-width NVT equilibration.
    observables.sampler:run(nvt_equilibrium_steps)

    local final_width = piston and piston.width or args.pore_width
    local analysis_domain, confinement_geometry = analysis_region(
        args,
        final_width,
        particle_walls,
        flat_walls,
        box
    )
    if args.sampling.structure > 0 then
        analysis_domain = slab_analysis.fit_region(
            args.slab.width,
            args.slab.axis,
            analysis_domain,
            box
        )
        local z_lower = analysis_domain.lowest_corner[3]
        confinement_geometry.analysis_lower = z_lower
        confinement_geometry.analysis_upper =
            z_lower + analysis_domain.length[3]
    end
    h5_metadata.write_confinement(file, confinement_geometry)

    -- Configure production observables.
    if args.sampling.structure > 0 then
        local write_structure_modes = structure_observables == "all"
        local wvs = write_structure_modes
            and wavevectors_mod.create(args, box)
            or nil
        local orientation_names = {
            parallel = "parallel",
            perpendicular = "normal"
        }

        if write_structure_modes and wvs then
            ssf_analysis.compute_ssf(args, group_fluid, wvs, file, {
                location = {"observables", "global", "structure"},
                orientation_names = orientation_names
            })
        end

        if not write_structure_modes or wvs then
            slab_analysis.slab_analysis(
                args.slab.width,
                args.slab.axis,
                particle_fluid,
                wvs,
                analysis_domain,
                file,
                box,
                args,
                {
                    density_modes = write_structure_modes,
                    ssf = write_structure_modes,
                    location = {
                        "observables",
                        "slabs",
                        args.slab.axis
                    },
                    geometry_location = {
                        "geometry",
                        "slabs",
                        args.slab.axis
                    },
                    orientation_names = orientation_names
                }
            )
        end
    end
    observables.sampler:run(production_steps)

    barostat_integrator:disconnect()

end

function define_args(parser)
    parser:add_argument("output,o", {
        type = "string",
        action = parser.action.substitute_date_time,
        default = data_dir
            .. "/production/npt"
            .. "_rho{fluid_density:g}"
            .. "_T{temperature:.2f}"
            .. "_%Y%m%d_%H%M%S",
        help = "basename of output files"
    })

    parser:add_argument("overwrite", {
        type = "boolean",
        default = false,
        help = "overwrite output file"
    })

    parser:add_argument("random-seed", {
        type = "integer",
        action = parser.action.random_seed,
        help = "seed for the random number generator"
    })

    parser:add_argument("particle_wall_input", {
        type = "string",
        default = data_dir .. "/wall/amorphous_symmetrized_W30_rho0.9.h5",
        help = "H5MD file containing particle walls"
    })

    parser:add_argument("pore_length", {
        type = "number",
        default = 10,
        help = "lateral pore length"
    })

    parser:add_argument("pore_width", {
        type = "number",
        required = true,
        help = "pore width in the z direction"
    })

    parser:add_argument("box_width", {
        type = "number",
        required = true,
        help = "fixed simulation box width in the z direction"
    })

    parser:add_argument("fluid_density", {
        type = "number",
        default = 0.01,
        help = "initial fluid number density"
    })

    parser:add_argument("temperature", {
        type = "number",
        default = 0.8,
        help = "target temperature"
    })

    parser:add_argument("rate", {
        type = "number",
        default = 1,
        help = "heat-bath collision rate"
    })

    parser:add_argument("cutoff", {
        type = "number",
        default = 3,
        help = "fluid and particle-wall LJ cutoff radius"
    })

    parser:add_argument("smoothing", {
        type = "number",
        default = 0.005,
        help = "potential cutoff smoothing parameter"
    })

    parser:add_argument("flat_cutoff", {
        type = "number",
        default = 2,
        help = "integrated flat-wall cutoff radius"
    })

    parser:add_argument("flat_offset", {
        type = "number",
        default = 1.5,
        help = "flat-wall padding beyond half the pore width"
    })

    parser:add_argument("equilibrium_time", {
        type = "number",
        default = 500,
        help = "initial fixed-wall NVT equilibration time"
    })

    parser:add_argument("barostat_time", {
        type = "number",
        default = 500,
        help = "NPAT barostat equilibration time"
    })

    parser:add_argument("nvt_equilibrium_time", {
        type = "number",
        default = 0,
        help = "fixed-width NVT equilibration time after the transition"
    })

    parser:add_argument("production_time", {
        type = "number",
        default = 500,
        help = "fixed-width NVT production time"
    })

    parser:add_argument("timestep", {
        type = "number",
        default = 0.001,
        help = "integration time step"
    })

    local barostat = parser:add_argument_group(
        "barostat",
        {help = "normal-pressure control for the pore width"}
    )

    barostat:add_argument("target_pressure", {
        type = "number",
        default = 0,
        help = "target normal pressure"
    })

    barostat:add_argument("mass", {
        type = "number",
        default = 1000,
        help = "piston mass"
    })

    barostat:add_argument("damping", {
        type = "number",
        default = 0,
        help = "linear damping coefficient for the piston velocity"
    })

    barostat:add_argument("max_width_change", {
        type = "number",
        default = 0.01,
        help = "maximum pore-width change per barostat update"
    })

    barostat:add_argument("interval", {
        type = "integer",
        default = 1,
        help = "barostat update interval in MD steps"
    })

    barostat:add_argument("transition_to_nvt", {
        type = "boolean",
        default = false,
        help = "switch to NVT near the mean equilibrated pore width"
    })

    barostat:add_argument("average_time", {
        type = "number",
        default = 10,
        help = "final barostat time window used to average the pore width"
    })

    barostat:add_argument("max_transition_time", {
        type = "number",
        default = 500,
        help = "maximum NPAT-to-NVT transition time"
    })

    barostat:add_argument("switch_tolerance", {
        type = "number",
        default = 0.02,
        help = "absolute pore-width tolerance for the NVT switch"
    })

    local weights = parser:add_argument_group(
        "weights",
        {help = "particle- and flat-wall interaction strengths"}
    )

    weights:add_argument("particle", {
        type = "number",
        default = 1,
        help = "particle-wall interaction weight"
    })

    weights:add_argument("flat", {
        type = "number",
        default = 0,
        help = "flat-wall interaction weight"
    })

    local sampling = parser:add_argument_group(
        "sampling",
        {help = "sampling intervals"}
    )

    sampling:add_argument("thermodynamics", {
        type = "integer",
        default = 100,
        help = "thermodynamic sampling interval"
    })
    sampling:add_argument("barostat", {
        type = "integer",
        default = 100,
        help = "barostat diagnostic sampling interval (0 disables output)"
    })
    sampling:add_argument("particles", {
        type = "integer",
        default = 0,
        help = "fluid position/velocity sampling interval (0 disables output)"
    })
    sampling:add_argument("structure", {
        type = "integer",
        default = 0,
        help = "density/structure sampling interval"
    })
    sampling:add_argument("structure_observables", {
        type = "string",
        default = "all",
        help = "structure output set: all or thermo"
    })
    local wavevector = parser:add_argument_group(
        "wavevector",
        {help = "wavevector shells in reciprocal space"}
    )
    observables.utility.wavevector.add_options(wavevector, {
        tolerance = 0.01,
        max_count = 7
    })
    observables.utility.semilog_grid.add_options(wavevector, {
        maximum = 5,
        decimation = 0
    })

    local slab = parser:add_argument_group(
        "slab",
        {help = "observables over slabs along one axis"}
    )
    slab:add_argument("width", {
        type = "number",
        default = 0.25,
        help = "slab width"
    })
    slab:add_argument("axis", {
        type = "string",
        default = "z",
        help = "axis normal to the slabs"
    })
end
