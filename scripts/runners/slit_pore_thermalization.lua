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
local arg_writer = require("sidecar")


function read_sample(args)
    local input = utility.assert_kwarg(args, "input")
    local location = utility.assert_kwarg(args, "location")
    local fields = args.fields      -- optional, see phase_space.reader()
    local step = -1                 -- use the last step stored in the file

    --local sidecar = arg_writer.write_args_from_output(args)
    --print("Writing args to " .. sidecar)

    -- open H5MD file for reading
    local file = readers.h5md({path = input})

    -- read edge vectors of simulation domain from file,
    -- only cuboid boxes are supported, so 'edges' must be a diagonal matrix
    local edges = mdsim.box.reader({file = file, location = location})

    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({edges = edges})

    local length = {edges[1][1], edges[2][2], edges[3][3]}

    -- construct a phase space reader and sample
    local reader, sample = observables.phase_space.reader({
        file = file, location = location, fields = fields
    })
    -- read phase space sample at given step in file
    reader:read_at_step(step)

    -- determine system parameters from phase space sample
    local nparticle = assert(sample.nparticle)
    local dimension = assert(sample.dimension)

    -- close H5MD file
    file:close()

    -- create system state
    local particle = mdsim.particle({dimension = dimension, particles = nparticle})

    -- use phase space "sampler" to set 'all' particle data (e.g., positions
    -- and velocities) by converting the sample obtained from the file
    local all_group = mdsim.particle_groups.all({particle = particle})
    observables.phase_space({box = box, group = all_group}):set(sample)

    return length, box, particle
end



function main(args)
    -- read obstacles particle data and box dimensions from input file
    local length, box, particle_obstacles = read_sample({
        input = args.input, location = {"particles", "obstacles"}, fields = {"position", "species"}
    })

    local dimension = #length
    local pore_length = {length[1], length[2], args.pore_width}
    local pore_volume = 100*100*15.95
    local pore_volume_func = function() return pore_volume end
    local fluid_density = args.fluid_density
    local wall_gap = 5
    local z_fraction = (args.pore_width - wall_gap) / length[3]
    local nparticle = math.ceil(pore_volume * fluid_density)

    -- select all obstacle particles to be used by the observables
    local group_obstacles = mdsim.particle_groups.all({particle = particle_obstacles, label = 'obstacles'})

    local particle_fluid = mdsim.particle({
        particles = nparticle,
        dimension = dimension,
        species = 2,
        label = 'fluid'
    })

    local group_fluid = mdsim.particle_groups.all({particle = particle_fluid})

    mdsim.positions.lattice({
        box = box,
        particle = particle_fluid,
        slab = {1, 1, z_fraction}
    }):set()

    -- set initial particle velocities
    local boltzmann = mdsim.velocities.boltzmann({
        particle = particle_fluid,
        temperature = args.temperature
    }):set()

    local particle_dict = {
        obstacles = particle_obstacles,
        fluid = particle_fluid
    }

    particle_fluid.data['species'] = utility.repeat_element(0, particle_fluid.nparticle)

    pair_potential.create_pair_forces(args.cutoff, args.smoothing, particle_dict, "fluid", "obstacles", box)

    -- add velocity-Verlet integrator with Boltzmann distribution
    local integrator = mdsim.integrators.verlet_nvt_boltzmann({
        box = box,
        particle = particle_fluid,
        timestep = args.timestep,
        temperature = args.temperature,
        rate = args.rate,
    })


    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = args.overwrite})

    local equilibrium_steps = math.ceil(args.equilibrium_time / args.timestep)

    -- write phase space trajectory to H5MD file
    local phase_space_obstacles = observables.phase_space({box = box, group = group_obstacles}):writer({file = file, fields = {"position", "image", "species", "velocity"}, every = equilibrium_steps})
    local phase_space_fluid     = observables.phase_space({box = box, group = group_fluid}):writer({file = file, fields = {"position", "image", "species",  "velocity"}, every = args.sampling.trajectory})

    --write thermodynamic variables to H5MD file
    observables.thermodynamics({box = box, group = group_fluid, volume = pore_volume_func}):writer({file = file,
    fields = { "potential_energy", "pressure", "temperature" , "center_of_mass_velocity", "stress_tensor", "virial"},
    every = args.sampling.trajectory})



    -- sample initial state
    observables.sampler:sample()
    -- estimate remaining runtime
    observables.runtime_estimate({steps = equilibrium_steps})
    -- run simulation
    observables.sampler:run(equilibrium_steps)


end

function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time
        , default ="/group/ag_compstatphys/data/tolga/simulation/slit_pore__equilibration_D{pore_width:g}_rho{fluid_density:g}_T{temperature:.2f}_%Y%m%d_%H%M%S"
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
    sampling:add_argument("trajectory", {type = "integer", default = 100000, help = "for trajectory"})
    sampling:add_argument("structure", {type = "integer", default = 100000, help = "for density modes, static structure factor"})
    sampling:add_argument("correlation", {type = "integer", default = 100, help = "for correlation functions"})

end
