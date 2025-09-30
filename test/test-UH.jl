# 小例子
a = 2.5
θ = 1.0      # 单位：小时
N = 24       # 单位线长度 24 小时
uh = gamma_unit_hydrograph(a, θ, N)

# 假设 RS 为 1 步冲击（冲击在第 5 小时）
RS = zeros(48)
RS[5] = 10.0
QS = route_with_uh(RS, uh)

println("uh sum = ", sum(uh))
println("QS[5:15] = ", QS[5:15])

dist = Gamma(a, θ)
x_range = 0:0.1:N
pdf_values = pdf.(dist, x_range)

 plot(x_range, pdf_values,
          title="Gamma Probability Density Function",
          xlabel="Time",
          ylabel="Probability Density",
          label="PDF",
          linewidth=2,
          color=:purple)
