--
-- custom module to set pair potentials
--

local mdsim = require("halmd.mdsim")
local M = {}


function M.create_pair_forces_init(cutoff, smoothing, particle, box)
    -- truncated Lennard-Jones potential, according to Lorentz-Bhertelot rules
    -- define pair potential
    local potential = mdsim.potentials.pair.lennard_jones({
        epsilon = {
            {1, 1} -- AA, AB
          , {1, 0} -- BA, BB
        }
      , sigma = {
            {1, 1} -- AA, AB
          , {1, 1} -- BA, BB
        }
    })
    -- apply smooth interaction cutoff
    if smoothing > 0 then
        potential = potential:truncate({"smooth_r4", cutoff = cutoff, h = smoothing})
    else
        potential = potential:truncate({cutoff = cutoff})
    end

    -- register computation of pair forces
    local force = mdsim.forces.pair({
        box = box, particle = particle, potential = potential
    })

    return force
end

function M.create_pair_forces(cutoff, smoothing, particle, label1, label2, box)
    -- truncated Lennard-Jones potential, according to Lorentz-Bhertelot rules
    local potential = mdsim.potentials.pair.lennard_jones({ 
        epsilon = { {1 , 1}, -- PP, PO
                    {1 , 0}  -- OP, OO
                    },
        sigma = { {1, 1},  -- PP, PO --sigma set to
                  {1, 1}   -- OP, OO
              }
        }):truncate({"smooth_r4", cutoff = cutoff, h = smoothing}) 
    --[[
    setting both σ and ϵ  to zero implies particles do not interact with each other at all

    σ=0: there is no finite distance at which the inter-particle potential becomes zero. Particles will interact even at infinitesimally small distances.
    ϵ=0: the potential well depth is zero, indicating that there are no attractive or repulsive forces between particles.
    --]]

     -- compute forces only based on the interaction between fluid-obstacles and fluid-fluid
    local fluid_obs = mdsim.forces.pair({
        box = box,
        particle = {particle[label1], particle[label2]},
        potential = potential
      , neighbour = { disable_sorting = true }
    })

    local fluid_fluid = mdsim.forces.pair({
        box = box,
        particle = particle[label1],
        potential = potential
    })

    return fluid_obs, fluid_fluid
end

return M