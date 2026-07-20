-- Simulate a homogeneous bulk state at a prescribed density.
local mdsim = halmd.mdsim
local numeric = halmd.numeric
local observables = halmd.observables
local writers = halmd.io.writers
local utility = halmd.utility

-- Load shared simulation helpers relative to this script.
package.path = utility.abspath("../../lib") .. "/?.lua;" .. package.path
local definitions = {lennard_jones = require("lennard_jones")}
local h5_metadata = require("h5_metadata")
local data_dir = utility.abspath("../../../data")


function main(args)
    local box_length = {30, 30, 30}
    local dimension = #box_length
    local density = args.density
    local volume = numeric.prod(box_length)
    local num_particles = math.ceil(volume * density)

    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({length = box_length})

    -- create system state
    local particle = mdsim.particle({
        dimension = dimension,
        particles = num_particles
    })

    -- set initial particle positions
    mdsim.positions.lattice({
        box = box,
        particle = particle
    }):set()

    -- set initial particle velocities
    local boltzmann = mdsim.velocities.boltzmann({
        particle = particle,
        temperature = args.temperature
    }):set()

    -- Define the Lennard-Jones pair force.
    definitions.lennard_jones.create_pair_force({
        box = box,
        particle = particle,
        cutoff = args.cutoff,
        smoothing = args.smoothing
    })

    -- convert integration time to number of steps
    local steps = math.ceil(args.time / args.timestep)

    -- H5MD file writer
    local file = writers.h5md({
        path = ("%s.h5"):format(args.output),
        overwrite = args.overwrite
    })

    h5_metadata.write_run(file, args)

    -- select all particles
    local group = mdsim.particle_groups.all({particle = particle})

    local grid = observables.utility.semilog_grid({
        start = 2 * math.pi / numeric.max(box.length),
        stop = args.wavevector.maximum,
        decimation = args.wavevector.decimation
    }).value

    local wavevector_parallel = observables.utility.wavevector({
        box = box,
        wavenumber = grid,
        tolerance = args.wavevector.tolerance,
        max_count = args.wavevector.max_count,
        filter = {1, 1, 0}
    })

    -- Compute global density modes
    local density_mode_parallel = observables.density_mode({
        group = group,
        wavevector = wavevector_parallel
    })

    -- Write global structure factors
    observables.ssf({
        density_mode = density_mode_parallel,
        norm = group.size
    }):writer({
        file = file,
        location = {"ssf", "parallel"},
        every = 1000
    })

    -- sample phase space
    local phase_space = observables.phase_space({
        box = box,
        group = group
    })

    -- write trajectory of particle groups to H5MD file
    local trajectory_interval = args.sampling.trajectory or steps
    if trajectory_interval > 0 then
        phase_space:writer({
            file = file,
            fields = {"position", "velocity"},
            every = trajectory_interval
        })
    end

    -- Sample macroscopic state variables.
    local state_interval = args.sampling.state_vars
    if state_interval > 0 then
        local msv = observables.thermodynamics({
            box = box,
            group = group
        })
        msv:writer({
            file = file,
            every = state_interval
        })
    end

    -- sample initial state
    observables.sampler:sample()

    -- add velocity-Verlet integrator with Boltzmann distribution
    local integrator = mdsim.integrators.verlet_nvt_boltzmann({
        box = box,
        particle = particle,
        timestep = args.timestep,
        temperature = args.temperature,
        rate = args.rate
    })

    -- estimate remaining runtime
    observables.runtime_estimate({steps = steps})

    -- run simulation
    observables.sampler:run(steps)
end

-- Parse command-line arguments.
function define_args(parser)
    parser:add_argument("output,o", {
        type = "string",
        action = parser.action.substitute_date_time,
        default = data_dir
            .. "/production/bulk"
            .. "_rho{density:g}"
            .. "_T{temperature:.2f}"
            .. "_rc{cutoff:g}"
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
        help = "seed for random number generator"
    })

    parser:add_argument("density", {
        type = "number",
        default = 0.75,
        help = "particle number density"
    })

    parser:add_argument("cutoff", {
        type = "number",
        default = 3.5,
        help = "potential cutoff radius"
    })

    parser:add_argument("smoothing", {
        type = "number",
        default = 0.005,
        help = "cutoff smoothing parameter"
    })

    parser:add_argument("temperature", {
        type = "number",
        default = 0.8,
        help = "initial system temperature"
    })

    parser:add_argument("rate", {
        type = "number",
        default = 2,
        help = "heat bath collision rate"
    })

    parser:add_argument("time", {
        type = "number",
        default = 100,
        help = "integration time"
    })

    parser:add_argument("timestep", {
        type = "number",
        default = 0.001,
        help = "integration time step"
    })

    local sampling = parser:add_argument_group(
        "sampling",
        {help = "sampling intervals (0: disabled)"}
    )

    sampling:add_argument("trajectory", {
        type = "integer",
        help = "for trajectory"
    })

    sampling:add_argument("state-vars", {
        type = "integer",
        default = 1000,
        help = "for state variables"
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
        maximum = 15,
        decimation = 0
    })
end
