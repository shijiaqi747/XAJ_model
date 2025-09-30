import Parameters: @with_kw
using Parameters
import FieldMetadata: @bounds, bounds, @units, units

abstract type AbstractModel{FT} end
abstract type AbstractState{FT} end
abstract type AbstractRouting{FT} end
abstract type AbstractSingleRouting{FT} <: AbstractRouting{FT} end
abstract type LinearReservoir{FT} <: AbstractSingleRouting{FT} end

struct RoutingVoid{FT} <: AbstractSingleRouting{FT} end

@bounds @units @with_kw mutable struct MultiRouting{FT} <: AbstractRouting{FT}
  RS::AbstractSingleRouting{FT} = RoutingVoid{FT}()
  RI::AbstractSingleRouting{FT} = RoutingVoid{FT}()
  RG::AbstractSingleRouting{FT} = RoutingVoid{FT}()
end


@bounds @units @with_kw mutable struct XAJ{FT} <: AbstractModel{FT}
  K::FT = 0.95 | (0.2, 1.5) | "-"       # 蒸发折算系数
  C::FT = 0.14 | (0.05, 0.20) | "-"     # 深层蒸发系数

  IM::FT = 0.1 | (0.01, 0.3) | "-"      # 不透水面积比例

  WUM::FT = 15.0 | (5.0, 20.0) | "mm"   # 上层土壤平均蓄水容量
  WLM::FT = 85.0 | (10.0, 90.0) | "mm"  # 下层
  WDM::FT = 20.0 | (10.0, 120.0) | "mm" # 深层

  SM::FT = 5.0 | (10.0, 60.0) | "mm"    # 自由水蓄水容量

  B::FT = 0.3 | (0.1, 0.6) | "-"        # 蓄水容量曲线的指数参数
  EX::FT = 0.7 | (0.5, 2.0) | "-"       # 自由水蓄水容量曲线指数

  KI::FT = 0.3 | (0.01, 0.7) | "-"      # 壤中流出流系数
  KG::FT = 0.3 | (0.01, 0.7) | "-"      # 地下水出流系数

  routing::MultiRouting{FT} = MultiRouting{FT}(
    RS=GammaUH{FT}(),                   # 伽马单位线
    RI=LinearReservoirSM{FT}(),         # 线性水库
    RG=LinearReservoirGM{FT}()          # 线性水库
  )
end



@with_kw mutable struct StateXAJ{FT} <: AbstractState{FT}
  ET::FT = 0.0
  EU::FT = 0.0
  EL::FT = 0.0
  ED::FT = 0.0
  PE::FT = 0.0    # 净雨, P - ET

  WU::FT = 0.0
  WL::FT = 2.2
  WD::FT = 20.0
  W::FT = WU + WL + WD

  FR::FT = 0.0
  R::FT = 0.0     # 透水界面径流
  R_IM::FT = 0.0  # 不透水界面径流
  RS::FT = 0.0
  RI::FT = 0.0
  RG::FT = 0.0
  S::FT = 0.0     # 自由水蓄量
end


# 输出结构体
@with_kw mutable struct OutputXAJ{FT<:Real}
    ntime::Int
    states::Vector{StateXAJ{FT}}
    R_sim::Vector{FT}
end

# 构造函数，根据 ntime 初始化
function Output(model::XAJ{FT}; ntime::Int) where {FT<:Real}
    OutputXAJ{FT}(
        ntime,
        [StateXAJ{FT}() for i in 1:ntime],  # 初始化每个时间步的状态
        zeros(FT, ntime)                     # 初始化 R_sim 为零
    )
end
