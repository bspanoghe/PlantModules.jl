# Plant compartment shapes used by hydraulic functional module #

"""
    Shape

An abstract type representing the physical shape of a compartment, used by the hydraulic functional module.
"""
abstract type Shape end

"""
    Sphere <: Shape

A spherical compartment shape, defined by a single dimension: the radius.
"""
struct Sphere <: Shape end

"""
    Cylinder <: Shape

A cylindrical compartment shape, defined by two dimensions: the radius and the length.
"""
struct Cylinder <: Shape end

"""
    Cuboid <: Shape

A cuboidal compartment shape, defined by three dimensions: the length, the width and the height.
"""
struct Cuboid <: Shape end

"""
    getdimensionality(s::Shape)

Return the number of dimensions that define shape `s`.
"""
getdimensionality(s::Shape) = error("Function not yet defined for shape $s")
getdimensionality(::Sphere) = 1
getdimensionality(::Cylinder) = 2
getdimensionality(::Cuboid) = 3

"""
    correctdimensionality(s::Shape, var)

Correct the dimensions of a variable `var` to match the number of dimensions of shape `s`.

If `var` is a scalar, return a vector of the length matching `s` with the value of `var` used for every dimension.
If `var` is an array, either return the same vector if the dimensions match or throw an error if they do not.
"""
function correctdimensionality(s::Shape, var)
    @info "A scalar value was found for a variable that needs to be defined for every dimensions of the shape (D, ϕ_D, ϵ_D)." *
    "This value will be used for every dimension." maxlog = 1
    return fill(var, getdimensionality(s))
end

function correctdimensionality(s::Shape, var::AbstractArray)
    if length(var) != getdimensionality(s)
        error("Volumes of shape $s must have $(getdimensionality(s)) dimension$(getdimensionality(s) > 1 ? "s" : "")" *
            " for all dimensional variables (by default: dimensions `D`, extensibility `ϕ_D`, elasticity `ϵ_D`).")
    end
    return var
end

"""
    volume(s::Shape, ::AbstractArray)

Calculate the volume of a given shape.
"""
volume(s::Shape, ::AbstractArray) = error("Function not yet defined for shape $s")
"""
    volume(::Sphere, D::AbstractArray)

Calculate the volume of a sphere.
"""
volume(::Sphere, D::AbstractArray) = 4/3 * pi * D[1]^3 # Write dimensions in the order: radius
"""
    volume(::Cylinder, D::AbstractArray)

Calculate the volume of a cylinder. Dimensions are assumed to be in the order: radius, length.
"""
volume(::Cylinder, D::AbstractArray) = D[1]^2 * pi * D[2] # Write dimensions in the order: radius - length
"""
    volume(::Cuboid, D::AbstractArray)

Calculate the volume of a cuboid.
"""
volume(::Cuboid, D::AbstractArray) = D[1] * D[2] * D[3] # Write dimensions in the order: length - width - height


"""
    cross_area(s::Shape, ::AbstractArray)

Calculate the cross-sectional area of a given shape.
"""
cross_area(s::Shape, ::AbstractArray) = error("Function not yet defined for shape $s")
"""
    cross_area(::Sphere, D::AbstractArray)

Calculate the cross-sectional area of a sphere, defined as the area of a circle with the same radius.
"""
cross_area(s::Sphere, D::AbstractArray) = D[1]^2 * pi # Write dimensions in the order: radius
"""
    cross_area(::Cylinder, D::AbstractArray)

Calculate the cross-sectional area of a cylinder, defined as its base area. Dimensions are assumed to be in the order: radius, length.
"""
cross_area(s::Cylinder, D::AbstractArray) = D[1]^2 * pi # Write dimensions in the order: radius - length
"""
    cross_area(::Cuboid, D::AbstractArray)

Calculate the cross-sectional area of a cuboid, defined as the product of its first two dimensions.
"""
cross_area(s::Cuboid, D::AbstractArray) = D[1] * D[2] # Write dimensions in the order: length - width - height

"""
    surface_area(s::Shape, ::AbstractArray)

Calculate the surface area of a given shape.
"""
surface_area(s::Shape, ::AbstractArray) = error("Function not yet defined for shape $s")
"""
    surface_area(::Sphere, D::AbstractArray)

Calculate the surface area of a sphere.
"""
surface_area(::Sphere, D::AbstractArray) = 4 * pi * D[1]^2
"""
    surface_area(::Cylinder, D::AbstractArray)

Calculate the surface area of a cylinder. Dimensions are assumed to be in the order: radius, length.
"""
surface_area(::Cylinder, D::AbstractArray) = 2 * (D[1]^2 * pi) + (2 * D[1] * pi) * D[2]
"""
    surface_area(::Cuboid, D::AbstractArray)

Calculate the surface area of a cuboid.
"""
surface_area(::Cuboid, D::AbstractArray) = 2 * (D[1]*D[2] + D[1]*D[3] + D[2]*D[3])