local M = {}

-- derive filename from args.output
local function derive_sidecar_path(output_path)
    local base = output_path:gsub("%.%w+$", "")
    return base .. "_args.lua"
end

function M.write_args_from_output(tbl)
    local out = tbl.output
    local fname = derive_sidecar_path(out)

    local f = assert(io.open(fname, "w"))
    for k, v in pairs(tbl) do
        local val
        if type(v) == "string" then
            val = string.format("%q", v)
        else
            val = tostring(v)
        end
        f:write(string.format("%s = %s,\n", k, val))
    end
    f:close()

    return fname
end

return M
