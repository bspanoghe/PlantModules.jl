"""
    alter_defaults(default_changes::NamedTuple; default_params::NamedTuple, default_u0s::NamedTuple)

Returns an altered version of default\\_params and default\\_u0s,
where parameter- or initial variable values that are defined in default_changes are set to this value.  
"""
function alter_defaults(default_changes::NamedTuple;
	default_params::NamedTuple = PlantModules.default_params, default_u0s::NamedTuple = PlantModules.default_u0s
	)

	altered_default_params = _alter_defaults(default_params, default_changes)
	altered_default_u0s = _alter_defaults(default_u0s, default_changes)

	return altered_default_params, altered_default_u0s
end

function _alter_defaults(nt, default_changes)
	NamedTuple{keys(nt)}(
		isempty(moduletuple) ? () : NamedTuple{keys(moduletuple)}(
			(
				get(default_changes, paramname, paramvalue)
				for (paramname, paramvalue) in zip(keys(moduletuple), values(moduletuple))
			)
		)
		for moduletuple in values(nt)
	)
end