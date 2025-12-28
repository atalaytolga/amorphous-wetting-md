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
local ssf_analysis = require("ssf_analysis")
local wavevectors_mod = require("wavevectors")

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


    local integrator = mdsim.integrators.verlet_nvt_boltzmann({
        box = box,
        particle = particles["fluid"],
        timestep = args.timestep,
        temperature = args.temperature,
        rate = args.rate,
    })
    local integration_steps = math.ceil(args.time / args.timestep)

    -- write phase space trajectory to H5MD file
    local phase_space_obstacles = observables.phase_space({box = box, group = groups["obstacles"]}):writer({file = file, fields = {"position", "image", "species"}, every = integration_steps})
    local phase_space_fluid     = observables.phase_space({box = box, group = groups["fluid"]}):writer({file = file, fields = {"position", "image", "species"}, every = args.sampling.trajectory})

    --write thermodynamic variables to H5MD file
    observables.thermodynamics({box = box, group = groups["fluid"], volume = pore_volume_func}):writer({file = file,
    fields = { "potential_energy", "pressure", "temperature" , "center_of_mass_velocity", "stress_tensor", "virial"},
    every = args.sampling.trajectory})

    local wvs = wavevectors_mod.create(args, box)
    if wvs then
        ssf_analysis.compute_ssf(args, groups["fluid"], wvs, file)

        if args.sampling.structure > 0 then
            slab_analysis.slab_analysis(
                  args.slab.width
                , args.slab.axis
                , particles["fluid"]
                , wvs
                , length
                , file
                , box
                , args
            )
        end
    end

    observables.sampler:sample()
    -- estimate remaining runtime
    observables.runtime_estimate({steps = integration_steps})
    -- run simulation
    observables.sampler:run(integration_steps)
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

    parser:add_argument("overwrite", {type = "boolean", default = false, help = "overwrite output file"})

    parser:add_argument('pore_width', {type = 'number', default = 15, help = 'pore width'})
    parser:add_argument('fluid_density', {type = 'number', default = 0.7, help = 'fluid density'})

    parser:add_argument("temperature", {type = "number", default = 1, help = "target temperature"})
    parser:add_argument("rate", {type = "number", default = 10, help = "heat bath collision rate"})
    parser:add_argument("time", {type = "number", default = 10, help = "simulation integration time"})
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
