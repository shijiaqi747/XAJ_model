using UnPack
# 三层蒸发模型
function cal_evaporation!(
    P::T,   # 降雨量，mm
    PET::T,  # 潜在蒸散发，mm
    state::StateXAJ{T},
    model::XAJ{T}
) where {T<:Real}
    (; WU, WL, WD) = state
    (; K, C, WLM) = model

    # 蒸发能力
    PET = PET * K

    EU = EL = ED = 0.0

    if WU + P >= PET
        EU = PET
    elseif WU + P < PET && WL >= C * WLM
        EU = WU + P
        EL = (PET - EU) * WL / WLM
    elseif WU + P < PET && C * (PET - (WU + P)) <= WL < C * WLM
        EU = WU + P
        EL = C * (PET - EU)
    else # WL < C(PET - EU)
        EU = WU + P
        EL = WL
        ED = C * (PET - EU) - EL
    end

    # 限制条件
    EU = min(EU, WU + P)
    EL = min(EL, WL)
    ED = min(ED, WD)
    
    # 总蒸发
    ET = EU + EL + ED
 
    PE = P - ET
    @pack! state = ET, EU, EL, ED, PE
end


# 产流计算
function cal_runoff!(
    state::StateXAJ{T},
    model::XAJ{T}
) where {T<:Real}

    (; IM, WUM, WLM, WDM, B) = model
    (; PE, W) = state

    # 初始化结果
    R = 0.0
    R_IM = 0.0
    FR = 0.0

    if PE > 0
        # 不透水区产流
        R_IM = IM * PE

        # 流域总储水量
        WM = WUM + WLM + WDM
        WMM = WM * (1 + B)

        # 当前平均含水量比例
        θ = clamp(W / WM, 0.0, 1.0 - eps(T))  # 控制 θ 在 [0,1) 内
        a = WMM * (1 - (1 - θ)^(1 / (1 + B)))

        # 产流计算
        if PE + a < WMM
            term = (1 - (PE + a) / WMM)^(B + 1)
            R = PE - (WM - W) + WM * (1 - term)
        else
            R = PE - (WM - W)
        end

        R = max(0.0, R)

        # 更新产流面积比例
        FR = 1 - (1 - W / WM)^B
    end

    return state
    @pack! state = R, R_IM, FR
end

#水储量更新
function update_W!(state::StateXAJ{T}, model::XAJ{T}) where {T<:Real}
    (; WUM, WLM, WDM) = model
    (; WU, WL, WD, PE, R, EU, EL, ED) = state

    # 上层水分更新,“有降水 → 上层加水但不超过容量；无降水 → 减少或保持非负”
    if PE > 0
        WU = min(WU + PE - R, WUM)
    else
        WU = max(WU + PE, 0.0)
    end

    # 中层水分更新,“有降水 → 剩余水量给中层，受容量限制；无降水 → 扣掉蒸发 EL 并保持非负”
    if PE > 0
        WL = WU + WL + WD + PE - R - WU - WD
    else
        WL = max(WL - EL, 0.0)
    end
    WL = clamp(WL, 0.0, WLM)

    # 深层水分更新,“有降水 → 超过上层和中层容量给下层；无降水 → 扣掉深层蒸发 ED 并保持非负”
    if PE > 0
        WD = max(WU + WL + WD + PE - R - WUM - WLM, 0.0)
    else
        WD = max(WD - ED, 0.0)
    end
    WD = min(WD, WDM)
    W = WU + WL + WD

    @pack! state = WU, WL, WD, W
end


# 分水源
function divide_runoff!(state::StateXAJ{T}, FR1::T, model::XAJ{T}) where {T<:Real}
    (; R, PE, FR, S, R_IM) = state
    (; SM, EX, KI, KG) = model

    # 初始化
    FR = FR == 0.0 ? T(0.1) : FR
    S = S == 0.0 ? 0.5 * SM : S

    RS = RI = RG = T(0)

    if PE > 0
        # 产流面积比例
        FR = clamp(R / PE, T(0), T(1))

        # 自由水储量换算
        S = S * FR1 / FR
        S = clamp(S, T(0), SM)

        # 最大自由水容量
        MS = SM * (1 + EX)
        AU = MS * (1 - clamp(1 - S / SM, 0.0, 1.0)^(1 / (1 + EX)))

        # 地表径流计算
        if PE + AU < MS
            RS = FR * (PE + SM * (1 - (PE + AU)/MS)^(EX+1) - SM + S)
            S = S + (R - RS) / FR
        else
            RS = FR * (PE + S - SM)
            S = SM
        end

        RS = clamp(RS, T(0), R)
    else
        # PE <= 0 时的分水源计算
        WW = clamp(R / (SM * (1 + EX)), 0.0, 1.0)
        FR = 1 - (1 - WW)^(1 / (1 + EX))
    end

    # 壤中流和地下水
    RI = KI * S * FR
    RG = KG * S * FR

    # 更新剩余自由水储量
    S = S * (1 - KI - KG)

    # 累加不透水区产流
    RS =RS + R_IM

    @pack! state = FR, RS, RI, RG, S
end
