# 伽马单位线结构体
@bounds @units @with_kw mutable struct GammaUH{FT} <: AbstractSingleRouting{FT}
    α::FT = 2.5 | (1.0, 5.0) | "-"       # 伽马分布形状参数
    θ::FT = 1.0 | (0.1, 3.0) | "hr"      # 伽马分布尺度参数
    τ_max::Int = 10 | (5, 50) | "t"      # 单位线长度 (时间步)
    UH::Vector{FT} = Vector{FT}(undef, 0) # 初始化为空，稍后生成
end

# 使用 gamma_unit_hydrograph 初始化单位线
function init_gammaUH!(obj::GammaUH{FT}) where {FT<:Real}
    obj.UH .= gamma_unit_hydrograph(obj.α, obj.θ, obj.τ_max)
end


@bounds @units @with_kw mutable struct LinearReservoirSM{FT} <: LinearReservoir{FT}
  K::FT = 0.5 | (0.1, 5.0) | "hr"   # 蓄水常数 
  dt::FT = 1.0 | (0.1, 5.0) | "hr"  # 时间步长 (可选)
end

@bounds @units @with_kw mutable struct LinearReservoirGM{FT} <: LinearReservoir{FT}
  K::FT = 0.5 | (0.1, 5.0) | "hr"   # 蓄水常数 
  dt::FT = 1.0 | (0.1, 5.0) | "hr"  # 时间步长 (可选)
end

function init_routing!(model::XAJ{FT}) where {FT}
    # RS: 伽马单位线
    RS_obj = GammaUH{FT}()
    init_gammaUH!(RS_obj)
    model.routing.RS = RS_obj

    # RI: 壤中流线性水库
    model.routing.RI = LinearReservoirSM{FT}()

    # RG: 地下水线性水库
    model.routing.RG = LinearReservoirGM{FT}()
end


function route_xaj(states::Vector{StateXAJ{FT}}, routing::MultiRouting{FT}) where {FT}
    # 从 states 中提取三条产流序列
    RS_series = [s.RS for s in states]
    RI_series = [s.RI for s in states]
    RG_series = [s.RG for s in states]

    # 把 routing 解构为局部变量，代码更简洁
    (; RS, RI, RG) = routing


    QS = route_with_uh(RS_series, RS.UH)                       # 表层：伽马单位线
    QI = linear_reservoir(RI_series; K = RI.K, dt = RI.dt)    # 壤中流：线性水库
    QG = linear_reservoir(RG_series; K = RG.K, dt = RG.dt)    # 地下水：线性水库

    Qtot = QS .+ QI .+ QG

    return (QS = QS, QI = QI, QG = QG, Q = Qtot)
end
