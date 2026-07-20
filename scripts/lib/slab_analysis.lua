-- Shared slab-analysis helpers.
local mdsim = require("halmd.mdsim")
local observables = halmd.observables
local numeric = halmd.numeric
local h5 = assert(libhalmd.h5)

local M = {}
local axis_map = {x = 1, y = 2, z = 3}


local function clone(values)
    local result = {}
    for index, value in ipairs(values) do
        result[index] = value
    end
    return result
end


local function append_location(location, ...)
    local result = clone(location)
    for _, name in ipairs({...}) do
        result[#result + 1] = name
    end
    return result
end


local function normalize_region(region)
    if region.lowest_corner and region.length then
        return {
            lowest_corner = clone(region.lowest_corner),
            length = clone(region.length)
        }
    end

    -- Backward compatibility for callers that pass a centred length vector.
    local lowest_corner = {}
    for index, value in ipairs(region) do
        lowest_corner[index] = -0.5 * value
    end
    return {
        lowest_corner = lowest_corner,
        length = clone(region)
    }
end


function M.fit_region(slab_width, axis, region_spec, box)
    local axis_index = axis_map[axis]
    if not axis_index then
        error("invalid slab axis: " .. tostring(axis))
    end
    if slab_width <= 0 then
        error("slab width must be positive")
    end

    local region = normalize_region(region_spec)
    local lower = region.lowest_corner[axis_index]
    local physical_upper = lower + region.length[axis_index]
    if physical_upper <= lower then
        error("slab analysis region must have positive extent")
    end

    local box_corner = box:lowest_corner()
    local box_lower = box_corner[axis_index]
    local box_upper = box_lower + box.length[axis_index]
    local tolerance = 1e-10 * math.max(1, box.length[axis_index])

    if lower < box_lower - tolerance
        or physical_upper > box_upper + tolerance
    then
        error("slab analysis region exceeds the simulation box")
    end

    local physical_width = physical_upper - lower
    local slab_count = math.max(
        1,
        math.ceil(physical_width / slab_width - 1e-12)
    )
    local analysis_width = slab_count * slab_width

    if analysis_width > box.length[axis_index] + tolerance then
        error("full-width slabs do not fit inside the simulation box")
    end

    local padding = analysis_width - physical_width
    lower = lower - 0.5 * padding
    local upper = physical_upper + 0.5 * padding
    if lower < box_lower - tolerance or upper > box_upper + tolerance then
        error("symmetrically padded slabs do not fit inside the simulation box")
    end

    region.lowest_corner[axis_index] = lower
    region.length[axis_index] = analysis_width

    return region, slab_count
end


function M.define_slabs(slab_width, axis, region_spec, particle, box)
    local axis_index = axis_map[axis]
    if not axis_index then
        error("invalid slab axis: " .. tostring(axis))
    end
    if slab_width <= 0 then
        error("slab width must be positive")
    end

    local region = normalize_region(region_spec)
    local lower = region.lowest_corner[axis_index]
    local upper = lower + region.length[axis_index]
    if upper <= lower then
        error("slab analysis region must have positive extent")
    end

    local slab_count = math.max(
        1,
        math.ceil((upper - lower) / slab_width - 1e-12)
    )
    local slabs = {}

    for index = 1, slab_count do
        local slab_lower = lower + (index - 1) * slab_width
        local slab_upper = slab_lower + slab_width
        local lowest_corner = clone(region.lowest_corner)
        local length = clone(region.length)
        lowest_corner[axis_index] = slab_lower
        length[axis_index] = slab_width

        local geometry = mdsim.geometries.cuboid({
            lowest_corner = lowest_corner,
            length = length
        })
        local group = mdsim.particle_groups.region({
            particle = particle,
            geometry = geometry,
            selection = "included",
            box = box,
            label = string.format("region-%s%04d", axis, index)
        })

        slabs[#slabs + 1] = {
            group = group,
            lower_edge = slab_lower,
            upper_edge = slab_upper,
            center = 0.5 * (slab_lower + slab_upper),
            volume = numeric.prod(length)
        }
    end

    return slabs
end


local function write_slab_geometry(file, location, axis, slab_width, slabs)
    local lower_edge = {}
    local upper_edge = {}
    local center = {}
    local volume = {}

    for index, slab in ipairs(slabs) do
        lower_edge[index] = slab.lower_edge
        upper_edge[index] = slab.upper_edge
        center[index] = slab.center
        volume[index] = slab.volume
    end

    local writer = file:writer({location = location, mode = "truncate"})
    local group = writer.group
    group:write_attribute("axis", h5.string(), axis)
    group:write_attribute("requested_width", h5.float(), slab_width)
    group:write_attribute("count", h5.int(), #slabs)
    group:write_attribute("lower", h5.float(), slabs[1].lower_edge)
    group:write_attribute("upper", h5.float(), slabs[#slabs].upper_edge)

    group:write_attribute("lower_edge", h5.float_array(), lower_edge)
    group:write_attribute("upper_edge", h5.float_array(), upper_edge)
    group:write_attribute("center", h5.float_array(), center)
    group:write_attribute("volume", h5.float_array(), volume)
end


function M.slab_analysis(
    slab_width,
    axis,
    particle,
    wavevectors,
    region,
    file,
    box,
    args,
    options
)
    local interval = args.sampling.structure
    options = options or {}
    local write_density_modes = options.density_modes ~= false
    local write_ssf = options.ssf ~= false

    if (write_density_modes or write_ssf) and not wavevectors then
        error("wavevectors are required for slab density modes or SSF")
    end

    local slabs = M.define_slabs(slab_width, axis, region, particle, box)
    local scoped_location = options.location
    local orientation_names = options.orientation_names or {
        parallel = "parallel",
        perpendicular = "perpendicular"
    }

    if options.geometry_location then
        write_slab_geometry(
            file,
            options.geometry_location,
            axis,
            slab_width,
            slabs
        )
    end

    for index, slab in ipairs(slabs) do
        local legacy_name = "slab_" .. index
        local slab_name = string.format("slab_%04d", index)
        local slab_location = scoped_location
            and append_location(scoped_location, slab_name)
            or nil

        local modes = {}
        if write_density_modes or write_ssf then
            modes.parallel = observables.density_mode({
                group = slab.group,
                wavevector = wavevectors.parallel
            })
            modes.perpendicular = observables.density_mode({
                group = slab.group,
                wavevector = wavevectors.perpendicular
            })
        end

        for _, orientation in ipairs({"parallel", "perpendicular"}) do
            local output_name = orientation_names[orientation]
            local density_mode = modes[orientation]

            if write_density_modes then
                density_mode:writer({
                    file = file,
                    location = slab_location
                        and append_location(
                            slab_location,
                            "structure",
                            output_name,
                            "density_mode"
                        )
                        or {"density_mode", legacy_name, orientation},
                    every = interval
                })
            end

            if write_ssf then
                observables.ssf({
                    density_mode = density_mode,
                    norm = slab.group.size
                }):writer({
                    file = file,
                    location = slab_location
                        and append_location(
                            slab_location,
                            "structure",
                            output_name,
                            "static_structure_factor"
                        )
                        or {"ssf", legacy_name, orientation},
                    every = interval
                })
            end
        end

        observables.thermodynamics({
            box = box,
            group = slab.group,
            volume = slab.volume
        }):writer({
            file = file,
            location = slab_location
                and append_location(slab_location, "thermodynamics")
                or {"thermo", legacy_name},
            fields = {
                "potential_energy",
                "pressure",
                "temperature",
                "stress_tensor",
                "virial"
            },
            every = interval
        })
    end
end


return M
