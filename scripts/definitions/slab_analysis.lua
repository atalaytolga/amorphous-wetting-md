local mdsim = require("halmd.mdsim")
local utility = require("halmd.utility")
local observables = halmd.observables
local numeric = halmd.numeric

local M = {}

local function clone3(v)
  return {v[1], v[2], v[3]}
end

function M.define_slabs(slit_width, axis, length, particle, box)
    local geometry = {}
    local group = {}

    local axis_map = {x = 1, y = 2, z = 3}
    local axis_idx = axis_map[axis]

    if not axis_idx then
        error("Invalid axis: " .. tostring(axis))
    end

    local total_length = length[axis_idx]

    local N = math.ceil(total_length / slit_width)

    local zmin = -length[axis_idx]/2

    local half = {-length[1]/2, -length[2]/2, -length[3]/2}

    -- Set region_length (same as box length, but slice thickness along axis)
    local region_length = {length[1], length[2], length[3]}
    region_length[axis_idx] = slit_width

     local slab_volume = function() return numeric.prod(region_length) end


    for i = 1, N do
        local lowest_corner = clone3(half)
        lowest_corner[axis_idx] = zmin + (i - 1) * slit_width
    
        -- Create geometry
        geometry[i] = mdsim.geometries.cuboid({
            lowest_corner = lowest_corner,
            length = region_length
        })

        -- Create particle group
        group[i] = mdsim.particle_groups.region({
            particle = particle,
            geometry = geometry[i],
            selection = "included",
            box = box,
            label = string.format("region-%s%d", axis, i)
        })
    end

    return N, group, slab_volume
end


function M.slab_analysis(slit_width, axis, particle, wavevectors, length, file, box, args)
    local interval = args.sampling.structure

    -- Define regions (slabs along defined axis)
    local Nslabs, slab_group, slab_volume = M.define_slabs(slit_width, axis, length, particle, box)

    -- Loop over slabs
    for i = 1, Nslabs do

        local density_mode_parallel = observables.density_mode({
            group = slab_group[i],
            wavevector = wavevectors["parallel"]
        })

        local density_mode_perpendicular = observables.density_mode({
            group = slab_group[i],
            wavevector = wavevectors["perpendicular"]
        })

        -- Write per-slab density modes
        density_mode_parallel:writer({
            file = file,
            location = {"density_mode", "slab_" .. i, "parallel"},
            every = interval
        })

        density_mode_perpendicular:writer({
            file = file,
            location = {"density_mode", "slab_" .. i, "perpendicular"},
            every = interval
        })

        -- Write per-slab structure factors
        observables.ssf({
            density_mode = density_mode_parallel,
            norm = slab_group[i].size
        }):writer({
            file = file,
            location = {"ssf", "slab_" .. i, "parallel"},
            every = interval
        })

        observables.ssf({
            density_mode = density_mode_perpendicular,
            norm = slab_group[i].size
        }):writer({
            file = file,
            location = {"ssf", "slab_" .. i, "perpendicular"},
            every = interval
        })

        -- Write per-slab thermodynamic observables
        observables.thermodynamics({
            box = box,
            group = slab_group[i],
            volume = slab_volume,
        }):writer({
            file = file,
            location = {"thermo", "slab_" .. i},
            fields = { "potential_energy", "pressure", "temperature" , "stress_tensor", "virial"},
            every = interval
        })
    end

end
return M