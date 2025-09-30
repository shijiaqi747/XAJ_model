

function run_XAJ(Prcp::AbstractVector{T}, PET::AbstractVector{T};
                 model::XAJ{T},
                 state::StateXAJ{T}=StateXAJ{T}(),
                 output::OutputXAJ{T}=Output(model; ntime=length(Prcp))) where {T<:Real}

    ntime = length(Prcp)

    for t in 1:ntime

        # 上一时刻产流面积比例 FR1
        FR1 = t == 1 ? zero(T) : output.states[t-1].FR
        _P = Prcp[t]
        _PET = PET[t] * model.K

        # 不透水区产流
        state.R_IM = _P * model.IM
        P = _P * (1 - model.IM)

        # 蒸发
        cal_evaporation!(P, _PET,state, model)

        # 产流
        cal_runoff!(state, model)

        # 更新水储量
        update_W!(state, model)

        # 分水源
        divide_runoff!(state, FR1, model)

        # 保存当前状态
        output.states[t] = deepcopy(state)

        # 更新 FR1
        FR1 = state.FR
    end

    # 汇流
    res = route_xaj(output.states, model.routing)
    output.R_sim = res.Q  # 总径流

    return output
end


function output_to_df(output::OutputXAJ)
    ntime = output.ntime
    states = output.states

    # 获取 StateXAJ 的所有字段名
    fnames = fieldnames(StateXAJ{eltype(states[1].ET)})

    # 为每个字段生成列
    columns = Dict()
    for fname in fnames
        columns[fname] = [getfield(s, fname) for s in states]
    end

    # 添加总径流列
    columns[:R_sim] = output.R_sim

    # 转为 DataFrame
    df = DataFrame(columns)

    return df
end