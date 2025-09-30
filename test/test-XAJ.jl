using HydroModels, Test
using RTableTools, DataFrames

d = fread("./data/dat_Table2-11.csv")
replace_missing!(d, 0.0)

P = d.P
PET = d.Eo

state = StateXAJ(; WU=0.0, WL=2.2, WD=20.0)

model = XAJ{Float64}(; K=0.95, C=0.14, B=0.3,
    WUM=15.0, WLM=85.0, WDM=20.0, IM=0.0)

    
r = run_XAJ(P, PET; model, state)
