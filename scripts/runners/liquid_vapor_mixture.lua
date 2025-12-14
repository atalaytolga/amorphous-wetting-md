-- grab modules
local log = halmd.io.log
local mdsim = halmd.mdsim
local numeric = halmd.numeric
local observables = halmd.observables
local readers = halmd.io.readers
local writers = halmd.io.writers
local utility = halmd.utility

-- search definition files in the top-level path relative to the simulation script
package.path = utility.abspath("../../definitions/?.lua;") .. package.path


function read_sample(args)
    -- open H5MD file for reading
    local file_liquid = readers.h5md({path = args.input.liquid})
    local file_vapor = readers.h5md({path = args.input.vapor})

    local particles = {}
    local groups = {}
    local edges

    local particles_list = {"liquid", "vapor"}

    -- construct a phase space reader and sample
    local reader_liquid, sample_liquid = observables.phase_space.reader({
        file = file_liquid
        , location = {"particles", "all"}
        , fields = {"position", "velocity"}
    })

    local reader_vapor, sample_vapor = observables.phase_space.reader({
        file = file_vapor
        , location = {"particles", "all"}
        , fields = {"position", "velocity"}
    })

    -- read phase space sample at last step in file
    log.info("number of liquid particles: %d", sample_liquid.nparticle)
    reader_liquid:read_at_step(-1)

    log.info("number of vapor particles: %d", sample_vapor.nparticle)
    reader_vapor:read_at_step(-1)

    -- read edge vectors of simulation domain from particle group
    edges = mdsim.box.reader({file = file_liquid, location = {"particles", "all"}})

    -- determine system parameters from phase space sample
    local nparticle_liquid = assert(sample_liquid.nparticle)
    local nparticle_vapor = assert(sample_vapor.nparticle)

    local dimension = assert(sample_liquid.dimension)
    local box = mdsim.box({edges = edges})

    -- create system state
    particles["liquid"] = mdsim.particle({dimension = dimension, particles = nparticle_liquid})
    groups["liquid"] = mdsim.particle_groups.all({particle = particles["liquid"], label = "liquid"})

    particles["vapor"] = mdsim.particle({dimension = dimension, particles = nparticle_vapor})
    groups["vapor"] = mdsim.particle_groups.all({particle = particles["vapor"], label = "vapor"})

    observables.phase_space({box = box, group = groups["liquid"]}):set(sample_liquid)
    observables.phase_space({box = box, group = groups["vapor"]}):set(sample_vapor)

    local length = {edges[1][1], edges[2][2], edges[3][3]}

    return particles, groups, length, box
end



function main(args)
    -- read obstacles particle data and box dimensions from input file
    local particles, groups, length, box = read_sample(args)

    -- H5MD file writer
    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = args.overwrite})


    local geometry = halmd.mdsim.geometries.cuboid({lowest_corner = {-10,-10,-15}, length = {20, 20, 30}})
    -- redefine new groups
    groups.liquid_region = halmd.mdsim.particle_groups.region({
        particle  = particles.liquid,
        label     = "liquid_region",
        geometry  = geometry,
        selection = "included",
    })

    groups.vapor_region = halmd.mdsim.particle_groups.region({
        particle  = particles.vapor,
        label     = "vapor_region",
        geometry  = geometry,
        selection = "excluded",

        
    })


    local particles_liquid_new = mdsim.particle({
        dimension = 3,
        particles = groups.liquid_region.size,
    })

    local particles_vapor_new = mdsim.particle({
        dimension = 3,
        particles = groups.vapor_region.size,
    })

    groups.liquid_region:to_particle({ particle = particles_liquid_new, label = "liquid" })
    groups.vapor_region:to_particle({ particle = particles_vapor_new, label = "vapor" })

    -- define “all” groups on the new particle sets
    local group_liquid_new = mdsim.particle_groups.all({
    particle = particles_liquid_new,
    label    = "liquid",
    })

    local group_vapor_new = mdsim.particle_groups.all({
    particle = particles_vapor_new,
    label    = "vapor",
    })

    particles_liquid_new.data['species'] = utility.repeat_element(1, particles_liquid_new.nparticle)
    particles_vapor_new.data['species'] = utility.repeat_element(0, particles_vapor_new.nparticle)


    local equilibrium_steps = math.ceil(args.equilibrium_time / args.timestep)

    -- write phase space trajectory to H5MD file
    local phase_space_liquid = observables.phase_space({box = box, group = group_liquid_new}):writer({file = file, fields = {"position", "image"}, every = args.sampling.trajectory})
    local phase_space_vapor     = observables.phase_space({box = box, group = group_vapor_new}):writer({file = file, fields = {"position", "image"}, every = args.sampling.trajectory})
    
    -- sample initial state
    observables.sampler:sample()

    
end


function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time
        , default ="/group/ag_compstatphys/data/tolga/simulation/liquid_vapor_mixture_NVT_T{temperature:.2f}_%Y%m%d_%H%M%S"
        , help = "basename of output files"
    })


    local input = parser:add_argument_group("input", {help = "sampling intervals (0: disabled)"})
    input:add_argument("vapor", {type = "string", required = true, action = function(args, key, value)
        readers.h5md.check(value)
        args[key] = value
    end, help = "H5MD input file"})
    input:add_argument("liquid", {type = "string", required = true, action = function(args, key, value)
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
