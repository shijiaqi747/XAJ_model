# 三层蒸发模型
function cal_evaporation!(   
    P::FT,   # 降雨量，mm
    PET::FT,  # 潜在蒸散发，mm
    state::StateXAJ{FT},
    model::XAJ{FT}
) where {FT}
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

   @pack! state = ET, EU, EL, ED, PE
end


# 产流计算
function cal_runoff!(
    state::StateXAJ{FT},
    model::XAJ{FT}
) where {FT}

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
        θ = clamp(W / WM, 0.0, 1.0 - eps(FT))  # 控制 θ 在 [0,1) 内
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
function update_W!(state::StateXAJ{FT}, model::XAJ{FT}) where {FT}
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
function divide_runoff!(state::StateXAJ{FT}, model::XAJ{FT}) where {FT}
    (; R, PE, FR, S) = state  
    (; SM, EX, KI, KG) = model

    # 初始化 FR 和 S
    FR = FR == 0.0 ? 0.1 : FR
    S  = S == 0.0 ? 0.5 * SM : S

    # 产流面积比例
    if PE > 0
        FR = R / PE
    end

    # 自由水储量换算（上一时段相关）
    SS = S
    if R > 0 && FR1 > 0
        SS = FR1 * S / FR
    end

    # 最大自由水容量
    MS = SM * (1 + EX)

    # 可用自由水容量
    AU = MS * (1 - (1 - (SS * FR1) / (FR * SM))^(1 / (1 + EX)))

    # 地表径流计算
    if PE + AU < MS
        RS = FR * (PE - SM + SS + SM * (1 - (PE + AU) / SM)^(1 + EX))
    else
        RS = FR * (PE + SS - SM)
    end
    RS = min(RS, R)  # 径流不超过总产流 R

    # 更新自由水储量
    S = min((SS * FR1) / FR + (R - RS) / FR, SM)

    # 壤中流和地下水
    RI = KI * S * FR
    RG = KG * S * FR

    # 更新剩余自由水储量
    S = S * (1 - KI - KG)

    # 保存当前 FR 作为下一时段 FR1
    @pack! state = FR, RS, RI, RG, S
end
