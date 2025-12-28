-- grab modules
local log = halmd.io.log
local mdsim = halmd.mdsim
local numeric = halmd.numeric
local observables = halmd.observables
local readers = halmd.io.readers
local writers = halmd.io.writers
local utility = halmd.utility

-- search definition files in the top-level path relative to the simulation script
package.path = utility.abspath("../definitions") .. "/?.lua;" .. package.path

local pair_potential = require("pair_potential")
local slab_analysis = require("slab_analysis")

function read_sample(args)
    -- open H5MD file for reading
    local file = readers.h5md({path = args.input})

    local particles = {}
    local groups = {}
    local edges

    local particles_list = {"fluid", "obstacles"}


    for _, label in ipairs(particles_list) do

        -- construct a phase space reader and sample
        local reader, sample = observables.phase_space.reader({
            file = file
          , location = {"particles", label}
          , fields = {"position", "species", "velocity"}
        })
        --particles[label] = sample
        -- read phase space sample at last step in file
        log.info("number of %s particles: %d", label, sample.nparticle)
        reader:read_at_step(-1)
        -- read edge vectors of simulation domain from particle group
        edges = mdsim.box.reader({file = file, location = {"particles", label}})

        local box = mdsim.box({edges = edges})

        -- determine system parameters from phase space sample
        local nparticle = assert(sample.nparticle)
        local dimension = assert(sample.dimension)

        -- create system state
        particles[label] = mdsim.particle({dimension = dimension, particles = nparticle})
        groups[label] = mdsim.particle_groups.all({particle = particles[label], label = label})

        observables.phase_space({box = box, group = groups[label]}):set(sample)

    end
    --local dimension = assert(particles[next(particles)].dimension)

    -- close H5MD file
    file:close()

    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({edges = edges})

    local length = {edges[1][1], edges[2][2], edges[3][3]}

    return particles, groups, length, box

end



function main(args)
    -- read obstacles particle data and box dimensions from input file
    local particles, groups, length, box = read_sample(args)

    local pore_length = {length[1], length[2], args.pore_width}
    local pore_volume = numeric.prod(pore_length)
    local pore_volume_func = function() return pore_volume end

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = args.overwrite})

    pair_potential.create_pair_forces(args.cutoff, args.smoothing, particles, "fluid", "obstacles", box)
    --ssf_analysis.compute_ssf(args, box, groups, file)


    local integrator = mdsim.integrators.verlet_nvt_boltzmann({
        box = box,
        particle = particles["fluid"],
        timestep = args.timestep,
        temperature = args.temperature,
        rate = args.rate,
    })

    --[[ add velocity-Verlet integrator with Boltzmann distribution
    -- add velocity-Verlet integrator (NVE)
    mdsim.integrators.verlet({particle = particles["fluid"], box = box, timestep = args.timestep})
    --]]
    local equilibrium_steps = math.ceil(args.equilibrium_time / args.timestep)

    -- write phase space trajectory to H5MD file
    local phase_space_obstacles = observables.phase_space({box = box, group = groups["obstacles"]}):writer({file = file, fields = {"position", "image", "species"}, every = equilibrium_steps})
    local phase_space_fluid     = observables.phase_space({box = box, group = groups["fluid"]}):writer({file = file, fields = {"position", "image", "species"}, every = args.sampling.trajectory})

    --write thermodynamic variables to H5MD file
    observables.thermodynamics({box = box, group = groups["fluid"], volume = pore_volume_func}):writer({file = file,
    fields = { "potential_energy", "pressure", "temperature" , "center_of_mass_velocity", "stress_tensor", "virial"},
    every = args.sampling.trajectory})


        -- set up wavevectors, compute density modes and static structure factor
    local interval = args.sampling.structure
    if interval > 0 or args.sampling.correlation > 0 then
        -- set up wavevector grid compatible with the periodic simulation box
        local grid = args.wavevector.wavenumbers
        if not grid then
            grid = observables.utility.semilog_grid({
                start = 2 * math.pi / numeric.max(box.length)
              , stop = args.wavevector.maximum
              , decimation = args.wavevector.decimation
            }).value
        end

        local wavevector_parallel = observables.utility.wavevector({
            box = box, wavenumber = grid
          , tolerance = args.wavevector.tolerance, max_count = args.wavevector.max_count
          , filter = {1,1,0}
        })

        local wavevector_perpendicular = observables.utility.wavevector({box = box,
        wavenumber = grid,
        dense = true,
        filter = {0,0,1}
        })

        local wavevectors = {
            parallel = wavevector_parallel,
            perpendicular = wavevector_perpendicular
        }


        -- Compute global density modes
        local density_mode_parallel_global = observables.density_mode({
            group = groups["fluid"],
            wavevector = wavevectors["parallel"]
        })

        local density_mode_perpendicular_global = observables.density_mode({
            group = groups["fluid"],
            wavevector = wavevectors["perpendicular"]
        })


        -- Write global density modes
        density_mode_parallel_global:writer({
        file = file,
        location = {"density_mode", "global_parallel", "bulk"},
        every = interval
        })


        density_mode_perpendicular_global:writer({
            file = file,
            location = {"density_mode", "global_perpendicular", "bulk"},
            every = interval
        })


        -- Write global structure factors
        observables.ssf({
            density_mode = density_mode_parallel_global,
            norm = groups["fluid"].size
        }):writer({
            file = file,
            location = {"ssf", "global_parallel", "bulk"},
            every = interval
        })
        observables.ssf({
            density_mode = density_mode_perpendicular_global,
            norm = groups["fluid"].size
        }):writer({
            file = file,
            location = {"ssf", "global_perpendicular", "bulk"},
            every = interval
        })
        slab_analysis.slab_analysis(args.slab.width, args.slab.axis, particles["fluid"], wavevectors, length, file, box, args)
    end


    -- sample initial state
    observables.sampler:sample()
    -- estimate remaining runtime
    observables.runtime_estimate({steps = equilibrium_steps})
    -- run simulation
    observables.sampler:run(equilibrium_steps)




end

function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time
        , default ="/group/ag_compstatphys/data/tolga/simulation/slit_pore__production_NVT_D{pore_width:g}_rho{fluid_density:g}_T{temperature:.2f}_%Y%m%d_%H%M%S"
        , help = "basename of output files"
    })

    parser:add_argument("input", {type = "string", required = true, action = function(args, key, value)
        readers.h5md.check(value)
        args[key] = value
    end, help = "H5MD input file"})

    parser:add_argument("overwrite", {type = "boolean", default = true, help = "overwrite output file"})
    parser:add_argument('pore_width', {type = 'number', default = 15, help = 'pore width'})
    parser:add_argument('wall_density', {type = 'number', default = 2.5, help = 'obstacle density'})
    parser:add_argument('fluid_density', {type = 'number', default = 0.7, help = 'fluid density'})
    parser:add_argument("temperature", {type = "number", default = 1, help = "target temperature"})
    parser:add_argument("rate", {type = "number", default = 10, help = "heat bath collision rate"})
    parser:add_argument("equilibrium_time", {type = "number", default = 10, help = "equilibruim integration time"})
    parser:add_argument("production_time", {type = "number", default = 10, help = "production integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.001, help = "integration time step"})
    parser:add_argument("smoothing", {type = "number", default = 0.005, help = "cutoff smoothing parameter"})
    parser:add_argument("cutoff", {type = "float32", default = 3, help = "potential cutoff radius"})


    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", default = 1000, help = "for trajectory"})
    sampling:add_argument("structure", {type = "integer", default = 1000, help = "for density modes, static structure factor"})
    sampling:add_argument("correlation", {type = "integer", default = 100, help = "for correlation functions"})

    local wavevector = parser:add_argument_group("wavevector", {help = "wavevector shells in reciprocal space"})
    observables.utility.wavevector.add_options(wavevector, {tolerance = 0.01, max_count = 7})
    observables.utility.semilog_grid.add_options(wavevector, {maximum = 15, decimation = 0})

    local slab = parser:add_argument_group("slab", {help = "ssf, density modes and observables over slabs along axis"})
    slab:add_argument("width", {type = "number", default = 0.5, help = "number of slabs"})
    slab:add_argument("axis", {type = "string", default = "z", help = "axis along slabs oriented"})

end
