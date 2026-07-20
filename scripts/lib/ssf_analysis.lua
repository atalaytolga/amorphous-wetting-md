local observables = halmd.observables

local M = {}


local function append_location(location, ...)
    local result = {}
    for index, name in ipairs(location) do
        result[index] = name
    end
    for _, name in ipairs({...}) do
        result[#result + 1] = name
    end
    return result
end


function M.compute_ssf(args, group, wavevectors, file, options)
    local interval = args.sampling.structure
    if interval <= 0 or not wavevectors then
        return
    end

    options = options or {}
    local location = options.location
    local orientation_names = options.orientation_names or {
        parallel = "parallel",
        perpendicular = "perpendicular"
    }

    local function compute_and_write(orientation, wavevector)
        if not wavevector then
            return
        end

        local output_name = orientation_names[orientation]
        local density_mode = observables.density_mode({
            group = group,
            wavevector = wavevector
        })

        density_mode:writer({
            file = file,
            location = location
                and append_location(location, output_name, "density_mode")
                or {"density_mode", "global_" .. orientation},
            every = interval
        })

        observables.ssf({
            density_mode = density_mode,
            norm = group.size
        }):writer({
            file = file,
            location = location
                and append_location(
                    location,
                    output_name,
                    "static_structure_factor"
                )
                or {"ssf", "global_" .. orientation},
            every = interval
        })
    end

    compute_and_write("parallel", wavevectors.parallel)
    compute_and_write("perpendicular", wavevectors.perpendicular)
end


return M
