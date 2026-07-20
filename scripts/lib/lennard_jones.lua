-- Shared Lennard-Jones force helpers.
-- Copyright © 2023 Felix Höfling
--
-- This file is part of HALMD.
--
-- HALMD is free software: you can redistribute it and/or modify
-- it under the terms of the GNU Lesser General Public License as
-- published by the Free Software Foundation, either version 3 of
-- the License, or (at your option) any later version.
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU Lesser General Public License for more details.
--
-- You should have received a copy of the GNU Lesser General
-- Public License along with this program.  If not, see
-- <http://www.gnu.org/licenses/>.
--

local mdsim = require("halmd.mdsim")
local utility = require("halmd.utility")

local M = {}

-- Create a possibly truncated single-species Lennard-Jones pair force.
-- args requires box and particle; cutoff and smoothing are optional.
function M.create_pair_force(args)
    local box = utility.assert_kwarg(args, "box")
    local particle = utility.assert_kwarg(args, "particle")
    local cutoff = args.cutoff
    local smoothing = args.smoothing or 0.005

    -- Define the Lennard-Jones pair potential.
    local potential = mdsim.potentials.pair.lennard_jones()
    -- Optionally apply an interaction cutoff.
    if cutoff then
        if smoothing > 0 then
            potential = potential:truncate({
                "smooth_r4",
                cutoff = cutoff,
                h = smoothing
            })
        else
            potential = potential:truncate({cutoff = cutoff})
        end
    end

    -- Register pair-force computation.
    local force = mdsim.forces.pair({
        box = box, particle = particle, potential = potential
    })

    return force, potential
end

return M
