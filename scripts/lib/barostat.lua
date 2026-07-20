local core = require("halmd.mdsim.core")
local mdsim = require("halmd.mdsim")
local observables = require("halmd.observables")

local h5 = assert(libhalmd.h5)

local M = {}
local Piston = {}
Piston.__index = Piston

local wall_force_labels = {
    top_particle = true,
    bottom_particle = true,
    top_flat = true,
    bottom_flat = true
}


local function assert_kwarg(args, name)
    local value = args[name]
    if value == nil then
        error(("missing argument '%s'"):format(name), 3)
    end
    return value
end


local function is_finite(value)
    return type(value) == "number"
        and value == value
        and value < math.huge
        and value > -math.huge
end


local function assert_finite(name, value)
    if not is_finite(value) then
        error(("%s must be finite: %s"):format(name, tostring(value)), 3)
    end
end


local function assert_positive(name, value)
    assert_finite(name, value)
    if value <= 0 then
        error(("%s must be positive"):format(name), 3)
    end
end

local function assert_pore_width(name, value, max_width)
    assert_positive(name, value)
    if value >= max_width then
        error(("%s (%g) reaches or exceeds maximum pore width (%g)")
            :format(name, value, max_width), 3)
    end
end


local function limit_width_change(self, old_width, width)
    local max_change = self.max_width_change
    if not max_change then
        return width, false
    end

    local delta = width - old_width
    if delta > max_change then
        return old_width + max_change, true
    end
    if delta < -max_change then
        return old_width - max_change, true
    end
    return width, false
end


local function limit_width_velocity(self, velocity, dt)
    local max_change = self.max_width_change
    if not max_change then
        return velocity
    end

    local max_velocity = max_change / dt
    if velocity > max_velocity then
        return max_velocity
    end
    if velocity < -max_velocity then
        return -max_velocity
    end
    return velocity
end


local function disconnect_all(connections)
    for index = #connections, 1, -1 do
        connections[index]:disconnect()
    end
end


local function new_force_state()
    return {
        modules = {},
        wall_force_z = {
            top_particle = 0,
            bottom_particle = 0,
            top_flat = 0,
            bottom_flat = 0
        },
        particle_generalized_force = 0,
        flat_generalized_force = 0,
        generalized_force = 0
    }
end


function Piston:apply_width(old_width, new_width)
    assert_pore_width("pore width", new_width, self.max_width)

    local half_shift = 0.5 * (new_width - old_width)

    if self.particle_walls then
        local top = self.particle_walls.top
        local bottom = self.particle_walls.bottom

        top.particle:shift_position(self.box, {0, 0, half_shift})
        bottom.particle:shift_position(self.box, {0, 0, -half_shift})
    end

    if self.flat_walls then
        local offset = self.flat_offset + 0.5 * new_width
        self.flat_walls.top:set_offset({offset})
        self.flat_walls.bottom:set_offset({offset})
    end
end


function Piston:generalized_force()
    local state = self.force_state
    local force = state.generalized_force

    if not is_finite(force) then
        error((
            "generalized force must be finite: total=%s, particle=%s, flat=%s"
        ):format(
            tostring(force),
            tostring(state.particle_generalized_force),
            tostring(state.flat_generalized_force)
        ), 2)
    end

    return force
end


function Piston:acceleration()
    local net_force = self:generalized_force()
        - self.target_pressure * self.area
        - self.damping * self.velocity

    return net_force / self.mass
end


-- First half-kick and width drift of velocity-Verlet.
function Piston:integrate(dt)
    local accel0 = self:acceleration()
    local velocity = self.velocity + 0.5 * dt * accel0
    local old_width = self.width
    local new_width, limited = limit_width_change(
        self, old_width, old_width + dt * velocity
    )

    if limited then
        velocity = (new_width - old_width) / dt
    end

    self.velocity = velocity
    self:apply_width(old_width, new_width)
    self.width = new_width
end


-- Second half-kick of velocity-Verlet after forces have been updated.
function Piston:finalize(dt)
    local accel1 = self:acceleration()
    self.velocity = limit_width_velocity(
        self,
        self.velocity + 0.5 * dt * accel1,
        dt
    )

    assert_finite("pore-width velocity", self.velocity)
end


function Piston:recover_wall_forces()
    local state = self.force_state
    local previous = 0

    -- Keep all four output channels present. Inactive wall types remain 0.
    for label in pairs(state.wall_force_z) do
        state.wall_force_z[label] = 0
    end

    -- A force module adds to the shared fluid-force array. Consecutive
    -- cumulative snapshots therefore isolate each module's contribution.
    for _, entry in ipairs(self.wall_forces) do
        local label = entry.label
        local module_state = state.modules[label]
        local current = module_state.cumulative_z

        if current == nil then
            error("missing force snapshot for " .. label, 2)
        end

        local contribution = current - previous
        previous = current

        if state.wall_force_z[label] ~= nil then
            -- Convert force on the fluid to reaction force on the wall.
            state.wall_force_z[label] = -contribution
        end
    end

    state.particle_generalized_force = 0.5 * (
        state.wall_force_z.top_particle
        - state.wall_force_z.bottom_particle
    )

    state.flat_generalized_force = 0.5 * (
        state.wall_force_z.top_flat
        - state.wall_force_z.bottom_flat
    )

    state.generalized_force =
        state.particle_generalized_force
        + state.flat_generalized_force
end


function Piston:connect_force_recovery()
    local state = self.force_state

    for _, entry in ipairs(self.wall_forces) do
        local label = entry.label
        local force = entry.force

        local module_state = {
            cumulative_z = nil
        }
        state.modules[label] = module_state

        self.force_connections[#self.force_connections + 1] =
            force:on_append_apply(function()
                module_state.cumulative_z =
                    self.thermodynamics:total_force()[3]
            end)
    end

    self.force_connections[#self.force_connections + 1] =
        self.particle:on_append_force(function()
            self:recover_wall_forces()
        end)

    -- Fixed-width equilibration may leave a valid cached force. Reapplying
    -- the current wall positions invalidates it without changing geometry.
    self:apply_width(self.width, self.width)
    self.thermodynamics:total_force()
    self:generalized_force()
end


function Piston:connect_core()
    local advance = false
    local dt = self.timestep * self.interval

    -- The fluid integrator was constructed first. Its append hook releases
    -- the old force before this hook shifts the walls.
    self.core_connections[#self.core_connections + 1] =
        core:on_append_integrate(function()
            self.step = self.step + 1
            advance = self.step % self.interval == 0

            if advance then
                self:integrate(dt)
            end
        end)

    self.core_connections[#self.core_connections + 1] =
        core:on_prepend_finalize(function()
            if not advance then
                return
            end

            -- Refresh once after the wall drift. The fluid finalizer reuses
            -- this force for its own second velocity-Verlet half-kick.
            self.thermodynamics:total_force()
            self:finalize(dt)

            if self.on_update then
                self.on_update(self, self.step)
            end
        end)
end


function Piston:connect_diagnostics()
    local args = self.diagnostics
    if not args or args.interval <= 0 then
        return
    end

    local location = {"observables", "barostat"}
    local metadata_writer = args.file:writer({
        location = location,
        mode = "truncate"
    })
    local metadata = metadata_writer.group

    metadata:write_attribute(
        "target_pressure",
        h5.float(),
        self.target_pressure
    )
    metadata:write_attribute("lateral_area", h5.float(), self.area)
    metadata:write_attribute("maximum_pore_width", h5.float(), self.max_width)
    metadata:write_attribute("piston_mass", h5.float(), self.mass)
    metadata:write_attribute("damping", h5.float(), self.damping)
    metadata:write_attribute(
        "barostat_interval",
        h5.int(),
        self.interval
    )
    metadata:write_attribute(
        "sampling_interval",
        h5.int(),
        args.interval
    )
    metadata:write_attribute(
        "particle_wall_weight",
        h5.float(),
        args.particle_wall_weight
    )
    metadata:write_attribute(
        "flat_wall_weight",
        h5.float(),
        args.flat_wall_weight
    )
    metadata:write_attribute(
        "controller_state_fields",
        h5.string_array(),
        {"pore_width", "pore_width_velocity", "normal_pressure"}
    )
    metadata:write_attribute(
        "pressure_component_fields",
        h5.string_array(),
        {
            "particle_normal_pressure",
            "flat_normal_pressure",
            "pressure_asymmetry"
        }
    )

    local function top_pressure()
        local state = self.force_state
        return (
            state.wall_force_z.top_particle
            + state.wall_force_z.top_flat
        ) / self.area
    end

    local function bottom_pressure()
        local state = self.force_state
        return -(
            state.wall_force_z.bottom_particle
            + state.wall_force_z.bottom_flat
        ) / self.area
    end

    -- HALMD's append writer cannot safely persist raw Lua scalar callbacks.
    -- Store the six diagnostics in two named vector fields of a one-particle
    -- carrier instead. The attributes above define each component exactly.
    local particle = mdsim.particle({
        particles = 1,
        dimension = 3,
        species = 1,
        label = "barostat_diagnostics"
    })
    local group = mdsim.particle_groups.all({
        particle = particle,
        label = "barostat_diagnostics"
    })

    local function sync()
        local state = self.force_state
        particle.data.position = {{
            self.width,
            self.velocity,
            state.generalized_force / self.area
        }}
        particle.data.velocity = {{
            state.particle_generalized_force / self.area,
            state.flat_generalized_force / self.area,
            top_pressure() - bottom_pressure()
        }}
    end

    sync()
    local writer = observables.phase_space({
        box = self.box,
        group = group
    }):writer({
        file = args.file,
        location = location,
        fields = {
            controller_state = "position",
            pressure_components = "velocity"
        },
        every = args.interval
    })
    self.diagnostic_sync_connection = writer:on_prepend_write(sync)
    self.diagnostic_particle = particle
    self.diagnostic_writer = writer
end


function Piston:connect()
    if self.state ~= "dormant" then
        error("barostat can only be connected once", 2)
    end

    self.state = "connecting"
    local success, message = pcall(function()
        self:connect_force_recovery()
        self:connect_core()
        self:connect_diagnostics()
    end)

    if not success then
        if self.diagnostic_sync_connection then
            self.diagnostic_sync_connection:disconnect()
            self.diagnostic_sync_connection = nil
        end
        if self.diagnostic_writer then
            self.diagnostic_writer:disconnect()
            self.diagnostic_writer = nil
        end
        self.diagnostic_particle = nil

        disconnect_all(self.core_connections)
        disconnect_all(self.force_connections)
        self.core_connections = {}
        self.force_connections = {}
        self.force_state = new_force_state()
        self.state = "failed"
        error(message, 2)
    end

    self.state = "connected"
end


function Piston:disconnect()
    if self.state ~= "connected" then
        error("barostat is not connected", 2)
    end

    if self.diagnostic_writer then
        self.diagnostic_writer:disconnect()
        self.diagnostic_writer = nil
    end
    if self.diagnostic_sync_connection then
        self.diagnostic_sync_connection:disconnect()
        self.diagnostic_sync_connection = nil
    end
    self.diagnostic_particle = nil

    disconnect_all(self.core_connections)
    disconnect_all(self.force_connections)
    self.core_connections = {}
    self.force_connections = {}
    self.state = "disconnected"
end


function M.create(args)
    local box = assert_kwarg(args, "box")
    local particle = assert_kwarg(args, "particle")
    local thermodynamics = assert_kwarg(args, "thermodynamics")
    local wall_forces = assert_kwarg(args, "wall_forces")
    local area = assert_kwarg(args, "area")
    local initial_width = assert_kwarg(args, "initial_width")
    local max_width = assert_kwarg(args, "max_width")
    local mass = assert_kwarg(args, "mass")
    local timestep = assert_kwarg(args, "timestep")

    local target_pressure = args.target_pressure or 0
    local damping = args.damping or 0
    local velocity = args.velocity or 0
    local interval = args.interval or 1
    local max_width_change = args.max_width_change
    local particle_walls = args.particle_walls
    local flat_walls = args.flat_walls
    local flat_offset = args.flat_offset or 0
    local on_update = args.on_update
    local diagnostics = args.diagnostics

    assert_positive("piston area", area)
    assert_positive("maximum pore width", max_width)
    assert_pore_width("initial pore width", initial_width, max_width)
    assert_positive("piston mass", mass)
    assert_positive("timestep", timestep)
    assert_positive("barostat interval", interval)
    assert_finite("target pressure", target_pressure)
    assert_finite("piston damping", damping)
    assert_finite("initial pore-width velocity", velocity)

    if interval % 1 ~= 0 then
        error("barostat interval must be an integer", 2)
    end
    if damping < 0 then
        error("piston damping must not be negative", 2)
    end
    if max_width_change ~= nil then
        assert_positive("maximum pore-width change", max_width_change)
    end
    if flat_walls then
        assert_finite("flat-wall offset", flat_offset)
    end
    if not particle_walls and not flat_walls then
        error("particle walls and/or flat walls are required", 2)
    end
    if type(wall_forces) ~= "table" or #wall_forces == 0 then
        error("wall_forces must be a non-empty table", 2)
    end

    local labels = {}
    for _, entry in ipairs(wall_forces) do
        if type(entry) ~= "table"
            or entry.label == nil
            or entry.force == nil
        then
            error("each wall force requires a label and force module", 2)
        end
        if not wall_force_labels[entry.label] then
            error("unknown wall force label: " .. tostring(entry.label), 2)
        end
        if labels[entry.label] then
            error("duplicate force label: " .. entry.label, 2)
        end
        labels[entry.label] = true
    end

    if on_update ~= nil and type(on_update) ~= "function" then
        error("on_update must be a function", 2)
    end

    if diagnostics then
        assert_kwarg(diagnostics, "file")
        local diagnostic_interval = assert_kwarg(diagnostics, "interval")
        local particle_wall_weight = assert_kwarg(
            diagnostics,
            "particle_wall_weight"
        )
        local flat_wall_weight = assert_kwarg(
            diagnostics,
            "flat_wall_weight"
        )
        assert_finite("diagnostic interval", diagnostic_interval)
        assert_finite("diagnostic particle-wall weight", particle_wall_weight)
        assert_finite("diagnostic flat-wall weight", flat_wall_weight)
        if diagnostic_interval < 0 or diagnostic_interval % 1 ~= 0 then
            error("diagnostic interval must be a non-negative integer", 2)
        end
        if particle_wall_weight < 0 or flat_wall_weight < 0 then
            error("diagnostic wall weights must not be negative", 2)
        end
    end

    return setmetatable({
        box = box,
        particle = particle,
        thermodynamics = thermodynamics,
        wall_forces = wall_forces,
        particle_walls = particle_walls,
        flat_walls = flat_walls,
        flat_offset = flat_offset,
        area = area,
        width = initial_width,
        max_width = max_width,
        velocity = velocity,
        target_pressure = target_pressure,
        mass = mass,
        damping = damping,
        max_width_change = max_width_change,
        timestep = timestep,
        interval = interval,
        on_update = on_update,
        diagnostics = diagnostics,
        force_state = new_force_state(),
        force_connections = {},
        core_connections = {},
        step = 0,
        state = "dormant"
    }, Piston)
end


return M
