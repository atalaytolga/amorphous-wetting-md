-- grab modules
local mdsim = halmd.mdsim
local numeric = halmd.numeric
local observables = halmd.observables
local writers = halmd.io.writers
local utility = halmd.utility

-- search definition files in the top-level path relative to the simulation script
package.path = utility.abspath("../../?.lua;") .. package.path
local definitions = { lennard_jones = require("definitions/lennard_jones") }

function main(args)
    local box_length = {40, 40, 160}
    local dimension = #box_length
    local vapor_density = 0.0185
    local volume = numeric.prod(box_length)
    local num_particles = volume * vapor_density

    -- create simulation domain with periodic boundary conditions
    local simulation_box = mdsim.box({length = box_length})

    -- create system state
    local vapor_particles = mdsim.particle({dimension = dimension, particles = num_particles})

    -- set initial particle positions
    mdsim.positions.lattice({box = simulation_box, particle = vapor_particles})
       :set()
    -- set initial particle velocities
    local boltzmann = mdsim.velocities.boltzmann({
        particle = vapor_particles
      , temperature = args.temperature
    }):set()

    -- define Lennard-Jones pair potential (with parameters ε=1 and σ=1 for a single species)
    -- and register computation of pair forces
    definitions.lennard_jones.create_pair_force({
        box = simulation_box, particle = vapor_particles, cutoff = args.cutoff, smoothing = args.smoothing
    })

    -- convert integration time to number of steps
    local steps = math.ceil(args.time / args.timestep)

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = args.overwrite})

    -- select all particles
    local all_group = mdsim.particle_groups.all({particle = vapor_particles})

    -- sample phase space
    local phase_space = observables.phase_space({box = simulation_box, group = all_group})
    -- write trajectory of particle groups to H5MD file
    local interval = args.sampling.trajectory or steps
    if interval > 0 then
        phase_space:writer({
            file = file, fields = {"position", "velocity"}, every = interval
        })
    end

        -- Sample macroscopic state variables.
    local msv
    local interval = args.sampling.state_vars
    if interval > 0 then
        msv = observables.thermodynamics({box = simulation_box, group = all_group})
        msv:writer({file = file, every = interval})
    end

    -- sample initial state
    observables.sampler:sample()

    -- add velocity-Verlet integrator with Boltzmann distribution
    local integrator = mdsim.integrators.verlet_nvt_boltzmann({
        box = simulation_box
      , particle = vapor_particles
      , timestep = args.timestep
      , temperature = args.temperature
      , rate = args.rate
    })

    -- estimate remaining runtime
    local runtime = observables.runtime_estimate({steps = steps})

    -- run simulation
    observables.sampler:run(steps)
end


--
-- Parse command-line arguments.
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time,
        default = "lennard_jones_equilibration_rc{cutoff:g}_rho{density:g}_T{temperature:.2f}_%Y%m%d_%H%M%S", help = "basename of output files"})
    parser:add_argument("overwrite", {type = "boolean", default = false, help = "overwrite output file"})

    parser:add_argument("random-seed", {type = "integer", action = parser.action.random_seed,
        help = "seed for random number generator"})

    parser:add_argument("density", {type = "number", default = 0.75, help = "particle number density"})

    parser:add_argument("cutoff", {type = "float32", default = 1.5, help = "potential cutoff radius"})
    parser:add_argument("smoothing", {type = "number", default = 0.005, help = "cutoff smoothing parameter"})
    parser:add_argument("masses", {type = "vector", dtype = "number", default = {1}, help = "particle masses"})
    parser:add_argument("temperature", {type = "number", default = 1.5, help = "initial system temperature"})
    parser:add_argument("rate", {type = "number", default = 2, help = "heat bath collision rate"})
    parser:add_argument("time", {type = "number", default = 100, help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.001, help = "integration time step"})

    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "for state variables"})
end