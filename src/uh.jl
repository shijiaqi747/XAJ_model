using SpecialFunctions: gamma
using Distributions

function gamma_unit_hydrograph(a::Real, θ::Real, N::Integer)
  dist = Gamma(a, θ)
  uh = zeros(eltype(a), N)
  prev = 0.0
  for i in 1:N
    cur = cdf(dist, i)    # P(T ≤ i)
    uh[i] = max(cur - prev, 0.0)  # 积分差，保证非负
    prev = cur
  end
  # 归一化以防数值误差
  s = sum(uh)
  s > 0 && (uh ./= s)

  return uh
end


function route_with_uh(RS::AbstractVector{T}, uh::AbstractVector{T}) where {T<:Real}
  Tlen = length(RS)
  N = length(uh)
  QS = zeros(Tlen)
  for t in 1:Tlen
    acc = zero(T)
    # sum i=1..min(N,t): RS[t-i+1] * uh[i]
    m = min(N, t)
    for i in 1:m
      acc += RS[t-i+1] * uh[i]
    end
    QS[t] = acc
  end
  return QS
end
