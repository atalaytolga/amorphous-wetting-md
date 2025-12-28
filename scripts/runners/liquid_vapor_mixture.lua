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
local definitions = { lennard_jones = require("lennard_jones") }
local wavevectors_mod = require("wavevectors")
local ssf_analysis = require("ssf_analysis")
local slab_analysis = require("slab_analysis")

function read_sample(args)
    -- open H5MD file for reading
    local file_sim = readers.h5md({path = args.input})

    -- construct a phase space reader and sample
    local reader, sample = observables.phase_space.reader({
        file = file_sim
        , location = {"particles", "all"}
        , fields = {"position", "velocity"}
    })


    -- read phase space sample at last step in file
    log.info("number of particles: %d", sample.nparticle)
    reader:read_at_step(-1)

    -- read edge vectors of simulation domain from particle group
    local edges = mdsim.box.reader({file = file_sim, location = {"particles", "all"}})
    local box = mdsim.box({edges = edges})
    local length = {edges[1][1], edges[2][2], edges[3][3]}

    -- determine system parameters from phase space sample
    local nparticle = assert(sample.nparticle)
    local dimension = assert(sample.dimension)

    -- create system state
    local particle = mdsim.particle({dimension = dimension, particles = nparticle})
    local group = mdsim.particle_groups.all({particle = particle, label = "all"})

    observables.phase_space({box = box, group = group}):set(sample)


    return particle, group, length, box
end



function main(args)
    -- read obstacles particle data and box dimensions from input file
    local particle, group, length, box = read_sample(args)

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = args.overwrite})

    -- define Lennard-Jones pair potential (with parameters ε=1 and σ=1 for a single species)
    -- and register computation of pair forces
    definitions.lennard_jones.create_pair_force({
        box = box, particle = particle, cutoff = args.cutoff, smoothing = args.smoothing
    })

    local steps = math.ceil(args.time / args.timestep)
    local interval = args.sampling.trajectory or steps

    -- write phase space of particle group to H5MD file
    if interval > 0 then
        observables.phase_space({box = box, group = group}):writer({
            file = file, fields = {"position", "velocity"}, every = interval
        })
    end

    --write thermodynamic variables to H5MD file
    observables.thermodynamics({box = box, group = group}):writer({file = file,
    fields = { "potential_energy", "pressure", "temperature" , "center_of_mass_velocity", "stress_tensor", "virial"},
    every = args.sampling.state_vars})

    -- Create Wavevectors
    local wvs = wavevectors_mod.create(args, box)

    if wvs then
        ssf_analysis.compute_ssf(args, group, wvs, file)

        if args.sampling.structure > 0 then
            slab_analysis.slab_analysis(
                args.slab.width,
                args.slab.axis,
                particle,
                wvs,
                length,
                file,
                box,
                args
            )
        end
    end
    -- Integrator: Nosé-Hoover
    local integrator = mdsim.integrators.verlet_nvt_hoover({
        box = box
      , particle = particle
      , timestep = args.timestep
      , temperature = args.temperature
      , resonance_frequency = args.resonance_frequency
    })


    observables.sampler:sample()
    observables.runtime_estimate({steps = steps})
    observables.sampler:run(steps)
end


function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time
        , default ="./data/raw/liquid_vapor_mixture_NVT_T{temperature:.2f}_%Y%m%d_%H%M%S"
        , help = "basename of output files"
    })


    parser:add_argument("input", {type = "string", required = true, action = function(args, key, value)
        readers.h5md.check(value)
        args[key] = value
    end, help = "H5MD input file"})

    parser:add_argument("overwrite", {type = "boolean", default = true, help = "overwrite output file"})
    parser:add_argument("temperature", {type = "number", default = 0.7, help = "target temperature"})
    parser:add_argument("resonance_frequency", {type = "number", default = 5, help = "Nosé-Hoover resonance frequency"})
    parser:add_argument("time", {type = "number", default = 100, help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.001, help = "integration time step"})
    parser:add_argument("smoothing", {type = "number", default = 0.005, help = "cutoff smoothing parameter"})
    parser:add_argument("cutoff", {type = "float32", default = 3, help = "potential cutoff radius"})


    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", default = 1000,  help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "for state variables"})
    sampling:add_argument("structure", {type = "integer", default = 1000, help = "for density modes and ssf"})
    sampling:add_argument("correlation", {type = "integer", default = 100, help = "for correlation functions"})

    local wavevector = parser:add_argument_group("wavevector", {help = "wavevector shells in reciprocal space"})
    observables.utility.wavevector.add_options(wavevector, {tolerance = 0.01, max_count = 7})
    observables.utility.semilog_grid.add_options(wavevector, {maximum = 15, decimation = 0})

    local slab = parser:add_argument_group("slab", {help = "ssf, density modes and observables over slabs along axis"})
    slab:add_argument("width", {type = "number", default = 0.5, help = "number of slabs"})
    slab:add_argument("axis", {type = "string", default = "z", help = "axis along slabs oriented"})

end
