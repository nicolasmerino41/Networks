num_species = 256
######################## COMPLEX RULES #############################
####################################################################
####################################################################
####################################################################
######################## DEFINING BASIC MYSTRUCTS256 METHODS ####################################
#################################################################################################
struct MyStructs256{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{num_species, T}
    b::T

    # Custom constructor for automatic sum calculation
    function MyStructs256(a::SVector{num_species, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor allowing manual setting of both `a` and `b`
    MyStructs256(a::SVector{num_species, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for MyStructs256
Base.zero(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{num_species, T}(fill(zero(T), num_species)), zero(T))
Base.zero(x::MyStructs256{T}) where {T <: AbstractFloat} = MyStructs256(SVector{num_species, T}(fill(zero(T), num_species)), zero(T))
Base.oneunit(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructsnum_species(fill(oneunit(T), num_species), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::MyStructs256, y::MyStructs256) = isless(x.b, y.b)
Base.isless(x::MyStructs256, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyStructs256, scalar::Real) = MyStructs256(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyStructs256, scalar::Real) = MyStructs256(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyStructs256, scalar::Real) = MyStructs256(x.a .- scalar, x.b - scalar * num_species)
Base.:+(x::MyStructs256, scalar::Real) = MyStructs256(x.a .+ scalar, x.b + scalar * num_species)
Base.:*(x::MyStructs256, y::AbstractVector) = MyStructs256(x.a .* SVector{num_species, Float64}(y), sum(x.a .* SVector{num_species, Float64}(y)))
Base.:/(x::MyStructs256, y::AbstractVector) = MyStructs256(x.a ./ SVector{num_species, Float64}(y), sum(x.a ./ SVector{num_species, Float64}(y)))
Base.:*(y::AbstractVector, x::MyStructs256) = MyStructs256(SVector{num_species, Float64}(y) .* x.a, sum(SVector{num_species, Float64}(y) .* x.a))
Base.:/(y::AbstractVector, x::MyStructs256) = MyStructs256(SVector{num_species, Float64}(y) ./ x.a, sum(SVector{num_species, Float64}(y) ./ x.a))


# Define what a NaN is for MyStructs256
Base.isnan(x::MyStructs256) = isnan(x.b) || any(isnan, x.a)

# Adding a method in the sum function for MyStructs256
function Base.sum(structs::MyStructs256...)
    # Sum the 'a' vectors
    summed_a = sum([s.a for s in structs])

    # Sum the 'b' values
    summed_b = sum([s.b for s in structs])

    # Create a new MyStructs256 instance with the summed results
    return MyStructs256(summed_a, summed_b)
end

# Adding a method to maximum
# Define maximum for MyStructs256
function Base.maximum(a::MyStructs256, b::MyStructs256)
    return MyStructs256(max.(a.a, b.a))
end

# Define maximum for MyStructs256 with a scalar
function Base.maximum(a::MyStructs256, b::AbstractFloat)
    return MyStructs256(max.(a.a, b))
end

# Define maximum for a scalar with MyStructs256
function Base.maximum(a::AbstractFloat, b::MyStructs256)
    return MyStructs256(max.(a, b.a))
end

# Define maximum for MyStructs256
function Base.maximum(a::MyStructs256)
    return maximum(a.a)
end

# Define maximum for a matrix of MyStructs256
function Base.maximum(a::Matrix{MyStructs256{AbstractFloat}})
    # Extract all `b` values from each MyStructs256 element in the matrix and find the maximum
    return maximum(map(x -> x.b, a))
end