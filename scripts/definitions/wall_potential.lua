--
-- custom module to create wall potentials
--
local mdsim = require("halmd.mdsim")
local M = {}

-- define the first wall of the slit of two planar walls made of integrated Lennard-Jones particles
function M.create_wall_force(box, particle, pore_diameter, smoothing, label1)
    
    local wall_potential = mdsim.potentials.external.planar_wall({
        offset = {(pore_diameter + 2) / 2, (pore_diameter + 2) / 2} ,
        surface_normal = {{0, 0, -1}, {0, 0, 1}},
        epsilon =1,
        sigma = 1,
        wetting = 0,
        cutoff = 2^(1/6),
        smoothing = smoothing,
        species = 2
    
    })
    local fluid_obs  = mdsim.forces.external({box = box, particle = particle[label1], potential = wall_potential })
    return fluid_obs
end

return M