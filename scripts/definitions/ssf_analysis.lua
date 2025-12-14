local numeric = halmd.numeric
local observables = halmd.observables
 
 local M = {}
 


 function M.compute_ssf(args, box, groups, file)
    -- set up wavevectors, compute density modes and static structure factor
        local interval = args.sampling.structure
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
                box = box, wavenumber = grid
            , tolerance = args.wavevector.tolerance, max_count = args.wavevector.max_count
            , filter = {1,1,0} 
            })

            local wavevector_perpendicular = observables.utility.wavevector({box = box,
            wavenumber = grid,
            dense = true,
            filter = {0,0,1}
            })

            local wavevectors = {
                parallel = wavevector_parallel,
                perpendicular = wavevector_perpendicular
            }


            -- Compute global density modes
            local density_mode_parallel_global = observables.density_mode({
                group = groups["fluid"],
                wavevector = wavevectors["parallel"]
            })

            local density_mode_perpendicular_global = observables.density_mode({
                group = groups["fluid"],
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


            -- Write global structure factors
            observables.ssf({
                density_mode = density_mode_parallel_global,
                norm = groups["fluid"].size
            }):writer({
                file = file,
                location = {"ssf", "global_parallel"},
                every = interval
            })
            observables.ssf({
                density_mode = density_mode_perpendicular_global,
                norm = groups["fluid"].size
            }):writer({
                file = file,
                location = {"ssf", "global_perpendicular"},
                every = interval
            })
            
        end
    end


return M

