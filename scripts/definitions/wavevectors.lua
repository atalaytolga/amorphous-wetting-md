-- definitions/wavevectors.lua

local numeric = halmd.numeric
local observables = halmd.observables

local M = {}

function M.create(args, box)
    if args.sampling.structure <= 0 and args.sampling.correlation <= 0 then
        return nil
    end

    local grid = args.wavevectors.wavenumbers
    if not grid then
        grid = observables.utility.semilog_grid({
              start = 2 * math.pi / numeric.max(box.length)
            , stop = args.wavevectors.maximum
            , decimation = args.wavevector.decimation
        }).value
    end

    local wavevector_parallel = observables.utility.wavevector({
          box = box
        , wavenumber = grid
        , tolerance = args.wavevector.tolerance
        , max_count = args.wavevector.max_count
        , filter = {1, 1, 0}
    })

    local wavevector_perpendicular = observables.utility.wavevector({
          box = box
        , wavenumber = grid
        , dense = true
        , filter = {0, 0, 1}
    })

    return {
        parallel = wavevector_parallel,
        perpendicular = wavevector_perpendicular
    }
end

return M


