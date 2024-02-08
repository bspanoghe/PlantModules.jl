"""
    PlantSystem

A plant system describing the plant's structural and functional modules and how they are connected.
This structure contains all information required to simulate the plant system using the `solve` function.

# Fields
- `model_defaults`: Model-wide default parameters.
- `module_defaults`: Module-specific default parameters.
- `module_coupling`: Coupling between functional and structural modules.
- `struct_connections`: One or more graphs specifying how the structural modules are connected, optionally alongside the connections between graphs.
- `func_connections`: Additional functional information about the connections between structural modules.
- `MTK_systems`: A vector of all structural modules represented as ModelingToolkit ODESystems.
- `MTK_connections`: A vector of all equations describing the functional connections as ModelingToolkit Equations.
- `MTK_u0`: A vector of all initial values.
"""
struct PlantSystem
	model_defaults
	module_defaults
	module_coupling
	struct_connections
	func_connections

    MTK_systems
    MTK_connections
    MTK_u0
    @doc """
        PlantSystem(; model_defaults, module_defaults, module_coupling, struct_connections, func_connections)

    The inner constructor. I don't know what else to say here.
    """
    function PlantSystem(; model_defaults, module_defaults, module_coupling, struct_connections, func_connections)
        MTK_systems = missing
        MTK_connections = missing
        MTK_u0 = missing

        new(model_defaults, module_defaults, module_coupling, struct_connections, func_connections,
            MTK_systems, MTK_connections, MTK_u0)
    end
end