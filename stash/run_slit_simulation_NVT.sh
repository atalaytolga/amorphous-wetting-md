halmd slit_simulation.lua -v \
    --input /group/ag_compstatphys/data/tolga/walls/crystalline_wall_200.h5 \
    --output /group/ag_compstatphys/data/tolga/simulation/sidecar_test \
    --pore_width 15 \
    --fluid_density 0.7 \
    --temperature 1.0 \
    --equilibrium_time 5 \
    --production_time 5 \
    --slab num_slabs 15 \
    --slab axis "z" 