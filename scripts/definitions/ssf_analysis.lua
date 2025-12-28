-- definitions/ssf_analysis.lua
local numeric = halmd.numeric
local observables = halmd.observables

 local M = {}

 --- Compute global density modes and static structure factors
 -- @param args: Parsed command-line arguments
 -- @param box: The simulation box
 -- @param group: The particle group to analyze
 -- @param file: The H5MD file writer
 -- @return table: Table containing the constructed wavevectors
 function M.compute_ssf(args, box, group, file)
    -- set up wavevectors, compute density modes and static structure factor
        local interval = args.sampling.structure
        local wavevectors = {}

        if interval > 0 or args.sampling.correlation > 0 then
            -- set up wavevector grid compatible with the periodic simulation box

            local grid = args.wavevector.wavenumbers
            if not grid then
                grid = observables.utility.semilog_grid({
                    start = 2 * math.pi / numeric.max(box.length)
                  , stop = args.wavevector.maximum
                  , decimation = args.wavevector.decimation
                }).value
            end

            local wavevector_parallel = observables.utility.wavevector({
                  box = box
                , wavenumber = grid
                , tolerance = args.wavevector.tolerance
                , max_count = args.wavevector.max_count
                , filter = {1,1,0}
            })

            local wavevector_perpendicular = observables.utility.wavevector({
                  box = box
                , wavenumber = grid
                , dense = true
                , filter = {0,0,1}
            })

            wavevectors = {
                parallel = wavevector_parallel,
                perpendicular = wavevector_perpendicular
            }


            -- Compute global density modes
            local density_mode_parallel_global = observables.density_mode({
                group = group,
                wavevector = wavevectors["parallel"]
            })

            local density_mode_perpendicular_global = observables.density_mode({
                group = group,
                wavevector = wavevectors["perpendicular"]
            })


            -- Write global density modes
            density_mode_parallel_global:writer({
            file = file,
            location = {"density_mode", "global_parallel"},
            every = interval
            })


            density_mode_perpendicular_global:writer({
                file = file,
                location = {"density_mode", "global_perpendicular"},
                every = interval
            })

            if interval > 0 then
                -- Write global structure factors
                observables.ssf({
                    density_mode = density_mode_parallel_global,
                    norm = group.size
                }):writer({
                    file = file,
                    location = {"ssf", "global_parallel"},
                    every = interval
                })
                observables.ssf({
                    density_mode = density_mode_perpendicular_global,
                    norm = group.size
                }):writer({
                    file = file,
                    location = {"ssf", "global_perpendicular"},
                    every = interval
                })

        end
    end

    return wavevectors
end

return M

