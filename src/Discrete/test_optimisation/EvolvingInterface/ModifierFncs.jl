"""
    Base.insert!(A::ElasticMatrix, col, xy)

Insert a two-element vector `xy` into column `col` of an `ElasticMatrix` `A` in-place.
"""
function Base.insert!(A::ElasticMatrix, col, xy)
    i1 = 2col - 1
    i2 = 2col
    data = A.data
    insert!(data, i1, xy[1])
    insert!(data, i2, xy[2])
    return A
end

"""
    Base.insert!(A::ElasticArray, col, val)

Insert `val` at position `col` in a 1D `ElasticArray` in-place.
"""
function Base.insert!(A::ElasticArray, col, val)
    data = A.data
    insert!(data, col, val)
    return A
end

"""
    Base.push!(A::ElasticArray, val)

Push `val` onto the end of a 1D `ElasticArray` in-place.
Works for any element type (String, Int8, Float64, etc.).
"""
function Base.push!(A::ElasticArray, val)
    data = A.data
    push!(data, val)
    return A
end

"""
    Base.deleteat!(A::ElasticMatrix, col)

Delete column `col` from an `ElasticMatrix` in-place.
"""
function Base.deleteat!(A::ElasticMatrix, col)
    data = A.data
    deleteat!(data, 2col - 1)
    deleteat!(data, 2col - 1)
    return A
end

"""
    Base.deleteat!(A::ElasticVector, col)

Delete element at position `col` from an `ElasticVector` in-place.
"""
function Base.deleteat!(A::ElasticVector, col)
    data = A.data
    deleteat!(data, col)
    return A
end
