# Plant compartment shapes used by hydraulic functional module #

"""
    Shape

An abstract type representing the physical shape of a compartment, used by the hydraulic functional module.
Subtypes should have an `ϕ_D` and a `ϵ_D` field, two vectors describing the dimensional extensibility and 
dimensional elastic modulus over the shape's defining dimensions respectively.
""" #! Pretty sure there's a better term for "defining dimensions"
abstract type Shape end

"""
    struct Sphere{T} <: Shape

A spherical compartment shape, defined by a single dimension: the radius.

# Fields
- `ϕ_D`: A one-dimensional vector containing the dimensional extensibility along the radius.
- `ϵ_D`: A one-dimensional vector containing the dimensional elastic modulus along the radius.
"""
struct Sphere{T} <: Shape
    ϕ_D::T
    ϵ_D::T
    function Sphere(ϕ_D::Vector, ϵ_D::Vector)
        length(ϕ_D) != 1 && error("An array of length $(length(ϕ_D)) was given for ϕ_D while length 1 was expected for shape Sphere.")
        length(ϵ_D) != 1 && error("An array of length $(length(ϵ_D)) was given for ϵ_D while length 1 was expected for shape Sphere.")

        return new{promote_type(typeof(ϕ_D), typeof(ϵ_D))}(promote(ϕ_D, ϵ_D)...)
    end
end

"""
    Sphere(ϕ_D::Number, ϵ_D::Number)

Convenience constructor for `Sphere`
"""
Sphere(ϕ_D::Number = 1, ϵ_D::Number = 1) = Sphere(fill(ϕ_D, 1), fill(ϵ_D, 1))

"""
    struct Cylinder{T} <: Shape

A cylindrical compartment shape, defined by two dimensions: the radius and the length.

# Fields
- `ϕ_D`: A two-dimensional vector containing the dimensional extensibility along the radius and the length.
- `ϵ_D`: A two-dimensional vector containing the dimensional elastic modulus along the radius and the length.
"""
struct Cylinder{T} <: Shape
    ϕ_D::T
    ϵ_D::T
    function Cylinder(ϕ_D::Vector, ϵ_D::Vector)
        length(ϕ_D) != 2 && error("An array of length $(length(ϕ_D)) was given for ϕ_D while length 2 was expected for shape Cylinder.")
        length(ϵ_D) != 2 && error("An array of length $(length(ϵ_D)) was given for ϵ_D while length 2 was expected for shape Cylinder.")

        return new{promote_type(typeof(ϕ_D), typeof(ϵ_D))}(promote(ϕ_D, ϵ_D)...)
    end
end

"""
    Cylinder(ϕ_D::Number, ϵ_D::Number)

Construct a `Cylinder` with equal dimensional elastic moduli and equal extensibilities over all dimensions. 
"""
Cylinder(ϕ_D::Number = 1, ϵ_D::Number = 1) = Cylinder(fill(ϕ_D, 2), fill(ϵ_D, 2))

"""
    struct Cuboid{T} <: Shape

A cuboidal compartment shape, defined by three dimensions: the length, the width and the height.

# Fields
- `ϕ_D`: A three-dimensional vector containing the dimensional extensibility along the length, the width and the height.
- `ϵ_D`: A three-dimensional vector containing the dimensional elastic modulus along the length, the width and the height.
"""
struct Cuboid{T} <: Shape
    ϕ_D::T
    ϵ_D::T
    function Cuboid(ϕ_D::Vector, ϵ_D::Vector)
        length(ϕ_D) != 3 && error("An array of length $(length(ϕ_D)) was given for ϕ_D while length 3 was expected for shape Cuboid.")
        length(ϵ_D) != 3 && error("An array of length $(length(ϵ_D)) was given for ϵ_D while length 3 was expected for shape Cuboid.")

        return new{promote_type(typeof(ϕ_D), typeof(ϵ_D))}(promote(ϕ_D, ϵ_D)...)
    end
end

"""
    Cuboid(ϕ_D::Number, ϕ_D::Number)

Construct a `Cuboid` with equal dimensional elastic moduli and equal extensibilities over all dimensions. 
"""
Cuboid(ϕ_D::Number = 1, ϵ_D::Number = 1) = Cuboid(fill(ϕ_D, 3), fill(ϵ_D, 3))

checkdimensions(s::Shape, D::AbstractArray) = error("Function not yet defined for shape $s")
checkdimensions(::Sphere, D::AbstractArray) = length(D) != 1 && error("Spheres must have one dimension")
checkdimensions(::Cylinder, D::AbstractArray) = length(D) != 2 && error("Cylinders must have two dimensions")
checkdimensions(::Cuboid, D::AbstractArray) = length(D) != 3 && error("Cuboids must have three dimension")

"""
    volume(s::Shape, ::AbstractArray)

Calculate the volume of a given shape.
"""
volume(s::Shape, ::AbstractArray) = error("Function not yet defined for shape $s")
"""
    volume(::Sphere, D::AbstractArray)

Calculate the volume of a sphere.
"""
volume(s::Sphere, D::AbstractArray) = (checkdimensions(s, D); 4/3 * pi * D[1]^3) # Write dimensions in the order: radius
"""
    volume(::Cylinder, D::AbstractArray)

Calculate the volume of a cylinder. Dimensions are assumed to be in the order: radius, length.
"""
volume(s::Cylinder, D::AbstractArray) = (checkdimensions(s, D); D[1]^2 * pi * D[2]) # Write dimensions in the order: radius - length
"""
    volume(::Cuboid, D::AbstractArray)

Calculate the volume of a cuboid.
"""
volume(s::Cuboid, D::AbstractArray) = (checkdimensions(s, D); D[1] * D[2] * D[3]) # Write dimensions in the order: length - width - height


"""
    cross_area(s::Shape, ::AbstractArray)

Calculate the cross-sectional area of a given shape.
"""
cross_area(s::Shape, ::AbstractArray) = error("Function not yet defined for shape $s")
"""
    cross_area(::Sphere, D::AbstractArray)

Calculate the cross-sectional area of a sphere, defined as the area of a circle with the same radius.
"""
cross_area(s::Sphere, D::AbstractArray) = (checkdimensions(s, D); D[1]^2 * pi) # Write dimensions in the order: radius
"""
    cross_area(::Cylinder, D::AbstractArray)

Calculate the cross-sectional area of a cylinder, defined as its base area. Dimensions are assumed to be in the order: radius, length.
"""
cross_area(s::Cylinder, D::AbstractArray) = (checkdimensions(s, D); D[1]^2 * pi) # Write dimensions in the order: radius - length
"""
    cross_area(::Cuboid, D::AbstractArray)

Calculate the cross-sectional area of a cuboid, defined as the product of its first two dimensions.
"""
cross_area(s::Cuboid, D::AbstractArray) = (checkdimensions(s, D); D[1] * D[2]) # Write dimensions in the order: length - width - height

"""
    surface_area(s::Shape, ::AbstractArray)

Calculate the surface area of a given shape.
"""
surface_area(s::Shape, ::AbstractArray) = error("Function not yet defined for shape $s")
"""
    surface_area(::Sphere, D::AbstractArray)

Calculate the surface area of a sphere.
"""
surface_area(::Sphere, D::AbstractArray) = (checkdimensions(s, D); 4 * pi * D[1]^2)
"""
    surface_area(::Cylinder, D::AbstractArray)

Calculate the surface area of a cylinder. Dimensions are assumed to be in the order: radius, length.
"""
surface_area(::Cylinder, D::AbstractArray) = (checkdimensions(s, D); 2 * (D[1]^2 * pi) + (2 * D[1] * pi) * D[2])
"""
    surface_area(::Cuboid, D::AbstractArray)

Calculate the surface area of a cuboid.
"""
surface_area(::Cuboid, D::AbstractArray) = 2 * (D[1]*D[2] + D[1]*D[3] + D[2]*D[3])