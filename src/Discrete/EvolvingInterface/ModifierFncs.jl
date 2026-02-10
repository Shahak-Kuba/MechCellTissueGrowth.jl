"""
    Base.insert!(A::ElasticMatrix, col, xy)

Insert a two-element vector `xy` into a specific column `col` of an `ElasticMatrix` `A`.

This function modifies the `ElasticMatrix` in place. It inserts the elements of `xy` into the underlying data of `A` at the positions corresponding to the specified column.

# Arguments
- `A::ElasticMatrix`: The `ElasticMatrix` object to be modified.
- `col`: The column index where the elements will be inserted.
- `xy`: A two-element vector containing the values to be inserted into the matrix.

# Returns
The modified `ElasticMatrix` `A`.
"""
function Base.insert!(A::ElasticMatrix, col, xy)
    i1 = 2col - 1
    i2 = 2col
    data = A.data
    insert!(data, i1, xy[1])
    insert!(data, i2, xy[2])
    return A
end

function Base.insert!(A::ElasticArray, col, val)
    data = A.data
    insert!(data, col, val)
    return A
end


function Base.push!(A::ElasticVector, val::String)
    data = A.data
    push!(data, val)
    return A
end

"""
    Base.deleteat!(A::ElasticMatrix, col)

Delete elements from a specific column `col` of an `ElasticMatrix` `A`.

This function modifies the `ElasticMatrix` in place. It removes the elements of the specified column from the underlying data of `A`.

# Arguments
- `A::ElasticMatrix`: The `ElasticMatrix` object to be modified.
- `col`: The column index from which the elements will be removed.

# Returns
The modified `ElasticMatrix` `A`.
"""
function Base.deleteat!(A::ElasticMatrix, col)
    data = A.data
    deleteat!(data,2col-1)
    deleteat!(data,2col-1)
    return A
end

function Base.deleteat!(A::ElasticVector, col)
    data = A.data
    deleteat!(data,col)
    return A
end