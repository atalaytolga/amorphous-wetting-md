-- grab modules
local mdsim = halmd.mdsim
local numeric = halmd.numeric
local observables = halmd.observables
local writers = halmd.io.writers
local utility = halmd.utility

-- search definition files in the top-level path relative to the simulation script
package.path = utility.abspath("../definitions/?.lua;") .. package.path

local definitions = {
    lennard_jones = require("lennard_jones")
}

function main(args)

    local length = {200, 200, (args.pore_width+15)}
    local dimension = #length
    local obstacle_density = args.obstacle_density
    local volume = numeric.prod(length)
    local nparticles = math.ceil(volume * obstacle_density)
    local pore_length = {length[1], length[2], args.pore_width}

    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({length = length})


    -- fill simulation volume with obstacle particles of density obstacle_density
    local particle_all = mdsim.particle({
        particles = nparticles,
        dimension = dimension,
        species = 1
    })


    mdsim.positions.lattice({
        box = box,
        particle = particle_all
    }):set()

    local group_all = mdsim.particle_groups.all({particle = particle_all})

    

    -- set initial particle velocities
    local boltzmann = mdsim.velocities.boltzmann({
        particle = particle_all,
        temperature = args.temperature
    }):set()


    -- define Lennard-Jones pair potential (with parameters ε=1 and σ=1 for a single species)
    -- and register computation of pair forces
    definitions.lennard_jones.create_pair_force({
        box = box, particle = particle_all, cutoff = args.cutoff, smoothing = args.smoothing
    })


     -- H5MD file writer
     local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = true})



    -- add velocity-Verlet integrator with Boltzmann distribution
    local integrator = mdsim.integrators.verlet_nvt_boltzmann({
        box = box,
        particle = particle_all,
        timestep = args.timestep,
        temperature = args.temperature,
        rate = args.rate
    })

    -- sample initial state
    observables.sampler:sample()
     -- steps
     local steps = math.ceil(args.time/args.timestep)
    -- estimate remaining runtime
    local runtime = observables.runtime_estimate({steps = steps})
    observables.sampler:run(steps)


    -- define pore geometry
    local pore_geometry = mdsim.geometries.cuboid({
        lowest_corner = {
            -pore_length[1]/2,
            -pore_length[2]/2,
            -args.pore_width/2
        },
        length = pore_length
    })

    local group_obstacles = mdsim.particle_groups.region({
            particle = particle_all, 
            geometry = pore_geometry, 
            selection = 'excluded', 
            label = 'obstacles',
            fluctuating = false,
            species = 1
    })
    particle_all.data['species'] = utility.repeat_element(1, particle_all.nparticle)



    local phase_space = observables.phase_space({box = box, group = group_obstacles}):writer({file = file, fields = {"position", "species"}, every = 1})
    observables.sampler:run(1)

end



function define_args(parser)
    parser:add_argument("output,o", {type = "string"
        , default ="/group/ag_compstatphys/data/tolga/walls/amorphous_walls_W{pore_width:g}_rho{obstacle_density:g}"
        , help = "basename of output files"
    })
    parser:add_argument("overwrite", {type = "boolean", default = false, help = "overwrite output file"})
    parser:add_argument('pore_width', {type = 'number', default = 15, help = 'pore width'})
    parser:add_argument("temperature", {type = "number", default = 1, help = "target temperature"})
    parser:add_argument("rate", {type = "number", default = 10, help = "heat bath collision rate"})
    parser:add_argument("time", {type = "number", default = 10 , help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.001, help = "integration time step"})
    parser:add_argument('obstacle_density', {type = 'number', default = 0.7, help = 'fluid density'})
    parser:add_argument("smoothing", {type = "number", default = 0.005, help = "cutoff smoothing parameter"})
    parser:add_argument("cutoff", {type = "float32", default = 3.5, help = "potential cutoff radius"})


    local wavevector = parser:add_argument_group("wavevector", {help = "wavevector shells in reciprocal space"})
    observables.utility.wavevector.add_options(wavevector, {tolerance = 0.01, max_count = 7})
    observables.utility.semilog_grid.add_options(wavevector, {maximum = 15, decimation = 0})
end