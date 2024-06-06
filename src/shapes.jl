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
struct Sphere{T}<:Shape
    ϵ_D::T
    ϕ_D::T
    function Sphere(; ϵ_D::Vector, ϕ_D::Vector)
        length(ϵ_D) != 1 && error("An array of length $(length(ϵ_D)) was given for ϵ_D while length 1 was expected for shape Sphere.")
        length(ϕ_D) != 1 && error("An array of length $(length(ϕ_D)) was given for ϕ_D while length 1 was expected for shape Sphere.")

        return new{typeof(ϵ_D)}(ϵ_D, ϕ_D)
    end
end

"""
    Cilinder

A cilindrical compartment shape, defined by two dimensions: the radius and the length.

# Fields
- `ϵ_D`: A two-dimensional vector containing the dimensional elastic modulus along the radius and the length.
- `ϕ_D`: A two-dimensional vector containing the dimensional extensibility along the radius and the length.
"""
struct Cilinder{T}<:Shape
    ϵ_D::T
    ϕ_D::T
    function Cilinder(; ϵ_D::Vector, ϕ_D::Vector)
        length(ϵ_D) != 2 && error("An array of length $(length(ϵ_D)) was given for ϵ_D while length 2 was expected for shape Cilinder.")
        length(ϕ_D) != 2 && error("An array of length $(length(ϕ_D)) was given for ϕ_D while length 2 was expected for shape Cilinder.")

        return new{typeof(ϵ_D)}(ϵ_D, ϕ_D)
    end
end

"""
    Cuboid

A cuboidal compartment shape, defined by three dimensions: the length, the width and the height.

# Fields
- `ϵ_D`: A three-dimensional vector containing the dimensional elastic modulus along the length, the width and the height.
- `ϕ_D`: A three-dimensional vector containing the dimensional extensibility along the length, the width and the height.
"""
struct Cuboid{T}<:Shape
    ϵ_D::T
    ϕ_D::T
    function Cuboid(; ϵ_D::Vector, ϕ_D::Vector)
        length(ϵ_D) != 3 && error("An array of length $(length(ϵ_D)) was given for ϵ_D while length 3 was expected for shape Cuboid.")
        length(ϕ_D) != 3 && error("An array of length $(length(ϕ_D)) was given for ϕ_D while length 3 was expected for shape Cuboid.")

        return new{typeof(ϵ_D)}(ϵ_D, ϕ_D)
    end
end

"""
    volume(s::Shape, ::AbstractArray)

Calculate the volume of a given shape.
"""
volume(s::Shape, ::AbstractArray) = error("Function `volume` is not defined for shape $s")
"""
    volume(::Sphere, D::AbstractArray)

Calculate the volume of a sphere.
"""
volume(::Sphere, D::AbstractArray) = 4/3 * pi * D[1]^3 # Write dimensions in the order: radius
"""
    volume(::Cilinder, D::AbstractArray)

Calculate the volume of a cilinder. Dimensions are assumed to be in the order radius, length.
"""
volume(::Cilinder, D::AbstractArray) = D[1]^2 * pi * D[2] # Write dimensions in the order: radius - length
"""
    volume(::Cuboid, D::AbstractArray)

Calculate the volume of a cuboid.
"""
volume(::Cuboid, D::AbstractArray) = D[1] * D[2] * D[3] # Write dimensions in the order: length - width - height

"""
    cross_area(s::Shape, ::AbstractArray)

Calculate the cross-sectional area of a given shape.
"""
cross_area(s::Shape, ::AbstractArray) = error("Function `cross_area` is not defined for shape $s")
"""
    cross_area(::Sphere, D::AbstractArray)

Calculate the cross-sectional area of a sphere, defined as the area of a circle with its radius.
"""
cross_area(::Sphere, D::AbstractArray) = D[1]^2 * pi # Write dimensions in the order: radius
"""
    cross_area(::Cilinder, D::AbstractArray)

Calculate the cross-sectional area of a cilinder, defined as its base area. Dimensions are assumed to be in the order radius, length.
"""
cross_area(::Cilinder, D::AbstractArray) = D[1]^2 * pi # Write dimensions in the order: radius - length
"""
    cross_area(::Cuboid, D::AbstractArray)

Calculate the cross-sectional area of a cuboid, defined as the product of its **last** two dimensions.
"""
cross_area(::Cuboid, D::AbstractArray) = D[2] * D[3] # Write dimensions in the order: length - width - height


"""
    surface_area(s::Shape, ::AbstractArray)

Calculate the surface area of a given shape.
"""
surface_area(s::Shape, ::AbstractArray) = error("Function `surface_area` is not defined for shape $s")
"""
    surface_area(::Sphere, D::AbstractArray)

Calculate the surface area of a sphere.
"""
surface_area(::Sphere, D::AbstractArray) = 4 * pi * D[1]^2
"""
    surface_area(::Cilinder, D::AbstractArray)

Calculate the surface area of a cilinder. Dimensions are assumed to be in the order radius, length.
"""
surface_area(::Cilinder, D::AbstractArray) = 2 * (D[1]^2 * pi) + (2 * D[1] * pi) * D[2]
"""
    surface_area(::Cuboid, D::AbstractArray)

Calculate the surface area of a cuboid.
"""
surface_area(::Cuboid, D::AbstractArray) = 2 * (D[1]*D[2] + D[1]*D[3] + D[2]*D[3])