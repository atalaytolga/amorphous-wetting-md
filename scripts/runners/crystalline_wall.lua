-- grab modules
local mdsim = halmd.mdsim
local numeric = halmd.numeric
local observables = halmd.observables
local writers = halmd.io.writers
local utility = halmd.utility

function main(args)

    local length = {100, 100, (args.pore_width +15)}
    local dimension = #length
    local pore_length = {length[1], length[2], args.pore_width}
    local obstacle_density = args.obstacle_density

    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({length = length})

    -- define pore geometry
    local pore_geometry = mdsim.geometries.cuboid({
        lowest_corner = {
            -pore_length[1]/2,
            -pore_length[2]/2,
            -args.pore_width/2
        },
        length = pore_length
    })

    -- fill simulation volume with obstacle particles of density obstacle_density
    local particle_obstacle = mdsim.particle({
        particles = math.ceil(numeric.prod(length) * obstacle_density),
        dimension = dimension,
        species = 1,
        label = 'obstacles'
    })

    mdsim.positions.lattice({
        box = box,
        particle = particle_obstacle
    }):set()

    local group_obstacles = mdsim.particle_groups.region({
            particle = particle_obstacle, 
            geometry = pore_geometry, 
            selection = 'excluded', 
            label = 'obstacles',
            fluctuating = false
    })

    particle_obstacle.data['species'] = utility.repeat_element(1, particle_obstacle.nparticle)



    local file = writers.h5md({path = ("%s.h5"):format(args.output), overwrite = args.overwrite})


    local grid = observables.utility.semilog_grid({
                start = 2 * math.pi / numeric.max(box.length)
              , stop = args.wavevector.maximum
              , decimation = args.wavevector.decimation
            }).value

    local wavevector_parallel = observables.utility.wavevector({
        box = box, wavenumber = grid
        , tolerance = args.wavevector.tolerance, max_count = args.wavevector.max_count
        , filter = {1,1,0} 
    })

    -- Compute global density modes
    local density_mode_parallel = observables.density_mode({
        group = group_obstacles,
        wavevector = wavevector_parallel
    })

    -- Write global structure factors
    observables.ssf({
        density_mode = density_mode_parallel,
        norm = group_obstacles.size
    }):writer({
        file = file,
        location = {"ssf", "parallel"},
        every = 1
    })







    observables.phase_space({box = box, group = group_obstacles}):writer({
        file = file, 
        fields = {"position", "species"},
        every = 1
    })

    observables.sampler:sample()

end


function define_args(parser)
    parser:add_argument("output,o", {type = "string"
        , default ="/group/ag_compstatphys/data/tolga/walls/crystalline_walls_100_W{pore_width:g}_rho{obstacle_density:g}_ssf"
        , help = "basename of output files"
    })
    parser:add_argument("overwrite", {type = "boolean", default = false, help = "overwrite output file"})
    parser:add_argument('pore_width', {type = 'number', default = 15, help = 'pore width'})
    parser:add_argument('obstacle_density', {type = 'number', default = 2.5, help = 'obstacle density'})
    
    local wavevector = parser:add_argument_group("wavevector", {help = "wavevector shells in reciprocal space"})
    observables.utility.wavevector.add_options(wavevector, {tolerance = 0.01, max_count = 7})
    observables.utility.semilog_grid.add_options(wavevector, {maximum = 100, decimation = 0})

end