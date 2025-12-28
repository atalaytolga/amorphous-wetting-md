-- definitions/ssf_analysis.lua
local observables = halmd.observables

 local M = {}

 --- Compute global density modes and static structure factors
 -- @param args: Parsed command-line arguments
 -- @param group: The particle group to analyze
 -- @param wavevectors: Table containing the constructed wavevectors
 -- @param file: The H5MD file writer
 function M.compute_ssf(args, group, wavevectors, file)
    local interval = args.sampling.structure

    if interval <= 0 or not wavevectors then return end

    local function compute_and_write(label, wv)
        if not wv then return end

            -- Compute global density modes
            local density_mode  = observables.density_mode({
                group = group,
                wavevector = wv
            })

            -- Write global density modes
            density_mode:writer({
            file = file,
            location = {"density_mode", "global_" .. label},
            every = interval
            })

            -- Write global structure factors
            observables.ssf({
                density_mode = density_mode,
                norm = group.size
            }):writer({
                file = file,
                location = {"ssf", "global_" .. label},
                every = interval
            })

    end

    compute_and_write("parallel", wavevectors.parallel)
    compute_and_write("perpendicular", wavevectors.perpendicular)
end

return M

