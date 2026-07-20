local h5 = assert(libhalmd.h5)

local M = {}

local function env(name, default)
    local value = os.getenv(name)
    if value == nil or value == "" then
        return default
    end
    return value
end

local function group(file, location)
    local writer = file:writer({
        location = location,
        mode = "truncate"
    })
    return assert(writer.group)
end

local function write_attr(g, name, value)
    if value == nil then
        return
    end

    if type(value) == "number" then
        g:write_attribute(name, h5.float(), value)
    elseif type(value) == "boolean" then
        g:write_attribute(name, h5.string(), tostring(value))
    elseif type(value) == "table" then
        g:write_attribute(name, h5.string_array(), value)
    else
        g:write_attribute(name, h5.string(), tostring(value))
    end
end

local function write_attrs(file, location, attrs)
    local g = group(file, location)
    for name, value in pairs(attrs) do
        write_attr(g, name, value)
    end
    return g
end

function M.write_run(file, args, opts)
    opts = opts or {}

    local run = write_attrs(file, {"parameters", "run"}, {
        schema_name = "amorphous-wetting-md",
        experiment_type = opts.experiment_type
            or env("AWMD_EXPERIMENT_TYPE", "unknown"),
        experiment_name = opts.experiment_name
            or env("AWMD_EXPERIMENT_NAME", "unknown"),
        simulation_script = opts.simulation_script
            or env("AWMD_SIMULATION_SCRIPT", "unknown"),
        output = args.output,
        created_at = env("AWMD_CREATED_AT", nil)
    })
    run:write_attribute(
        "schema_version",
        h5.int(),
        opts.schema_version or 1
    )

    write_attrs(file, {"parameters", "provenance"}, {
        repo_commit = env("AWMD_GIT_COMMIT", nil),
        repo_branch = env("AWMD_GIT_BRANCH", nil),
        repo_dirty = env("AWMD_GIT_DIRTY", nil),
        backend = env("AWMD_BACKEND", nil),
        docker_image = env("AWMD_DOCKER_IMAGE", nil)
    })
end

function M.write_confinement(file, geometry)
    local g = group(file, {"geometry", "confinement"})

    g:write_attribute("axis", h5.string(), geometry.axis)
    g:write_attribute("pore_width", h5.float(), geometry.pore_width)
    g:write_attribute(
        "analysis_lower",
        h5.float(),
        geometry.analysis_lower
    )
    g:write_attribute(
        "analysis_upper",
        h5.float(),
        geometry.analysis_upper
    )
    g:write_attribute(
        "analysis_width",
        h5.float(),
        geometry.analysis_upper - geometry.analysis_lower
    )
    g:write_attribute(
        "top_wall_thickness",
        h5.float(),
        geometry.top_wall_thickness
    )
    g:write_attribute(
        "bottom_wall_thickness",
        h5.float(),
        geometry.bottom_wall_thickness
    )
    g:write_attribute(
        "flat_wall_padding",
        h5.float(),
        geometry.flat_wall_padding
    )
    g:write_attribute(
        "particle_walls_enabled",
        h5.int(),
        geometry.particle_walls_enabled and 1 or 0
    )
    g:write_attribute(
        "flat_walls_enabled",
        h5.int(),
        geometry.flat_walls_enabled and 1 or 0
    )
end

return M
