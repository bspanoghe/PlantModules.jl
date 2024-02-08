# Plant compartment shapes used by hydraulic functional module #

"""
    Shape

An abstract type representing the physical shape of a compartment, used by the hydraulic functional module.
Subtypes should have an `ϵ_D` and a `ϕ_D` field, two vectors describing the dimensional elastic modulus and 
dimensional extensibility over the shape's defining dimensions respectively.
""" #! Pretty sure there's a better term for "defining dimensions"
abstract type Shape end

"""
    Sphere

A spherical compartment shape, defined by a single dimension: the radius.

# Fields
- `ϵ_D`: A one-dimensional vector containing the dimensional elastic modulus along the radius.
- `ϕ_D`: A one-dimensional vector containing the dimensional extensibility along the radius.
"""
struct Sphere<:Shape
    ϵ_D::Vector
    ϕ_D::Vector
    function Sphere(; ϵ_D::Vector, ϕ_D::Vector)
        length(ϵ_D) != 1 && error("An array of length $(length(ϵ_D)) was given for ϵ_D while length 1 was expected.")
        length(ϕ_D) != 1 && error("An array of length $(length(ϕ_D)) was given for ϕ_D while length 1 was expected.")

        return new(ϵ_D, ϕ_D)
    end
end

"""
    Cilinder

A cilindrical compartment shape, defined by two dimensions: the radius and the length.

# Fields
- `ϵ_D`: A two-dimensional vector containing the dimensional elastic modulus along the radius and the length.
- `ϕ_D`: A two-dimensional vector containing the dimensional extensibility along the radius and the length.
"""
struct Cilinder<:Shape
    ϵ_D::Vector
    ϕ_D::Vector
    function Cilinder(; ϵ_D::Vector, ϕ_D::Vector)
        length(ϵ_D) != 2 && error("An array of length $(length(ϵ_D)) was given for ϵ_D while length 2 was expected.")
        length(ϕ_D) != 2 && error("An array of length $(length(ϕ_D)) was given for ϕ_D while length 2 was expected.")

        return new(ϵ_D, ϕ_D)
    end
end

"""
    Cuboid

A cuboidal compartment shape, defined by three dimensions: the length, the width and the height.

# Fields
- `ϵ_D`: A three-dimensional vector containing the dimensional elastic modulus along the length, the width and the height.
- `ϕ_D`: A three-dimensional vector containing the dimensional extensibility along the length, the width and the height.
"""
struct Cuboid<:Shape
    ϵ_D::Vector
    ϕ_D::Vector
    function Cuboid(; ϵ_D::Vector, ϕ_D::Vector)
        length(ϵ_D) != 3 && error("An array of length $(length(ϵ_D)) was given for ϵ_D while length 3 was expected.")
        length(ϕ_D) != 3 && error("An array of length $(length(ϕ_D)) was given for ϕ_D while length 3 was expected.")

        return new(ϵ_D, ϕ_D)
    end
end

volume(s::Shape, D::AbstractArray) = error("Function volume is not defined for shape $s")
volume(::Sphere, D::AbstractArray) = 4/3 * pi * D[1]^3 # Write dimensions in the order: radius
volume(::Cilinder, D::AbstractArray) = D[1]^2 * pi * D[2] # Write dimensions in the order: radius - length
volume(::Cuboid, D::AbstractArray) = D[1] * D[2] * D[3] # Write dimensions in the order: length - width - height
