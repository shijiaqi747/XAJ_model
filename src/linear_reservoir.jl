"""
Inputs:
- I :: Vector{T}   入流序列
- K :: Real        水库常数
- Δt :: Real       时间步长 


Return:
- Q :: Vector{T}   出流序列
"""
function linear_reservoir(I::AbstractVector{T}; K::Real, dt::Real=1.0) where {T<:Real}
    ntime = length(I)
    Q = zeros(T, ntime)

    Cs = (2K - dt) / (2K + dt)
    Q[1] = I[1]

    @inbounds for t in 2:ntime
        Q[t] = Cs * Q[t-1] + (1 - Cs) * (I[t-1] + I[t]) / 2
    end

    return Q
end
