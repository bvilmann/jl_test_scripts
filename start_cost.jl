using JuMP, Cbc, CPLEX              # Optimization and modelling, Gurobi ?
using Plots, LaTeXStrings       # Plotting

# ========================= DATA =========================
# Power capacity
P_C  = [50 200; 25 250; 75 300; 100 400; 125 500; 150 600; 175 700; 200 800; 225 900; 100 1000;]
# Power demand
P_D = [3268.78 3158.69 3036.5 2946.15 2891.74 2909.41 3014.69 3104.31 3243.79 3547.45 3808.79 4134.67 4157.68 4048.36 4184.9 4052.32 4229.83 4510.66 4403.63 4158.34 3916.65 3669.91 3421.24 3164.72 ]*1
P_D = [P_D P_D P_D ]                            # Create 3 cycles of power demand
T = length(P_D)                                 # Number of time steps
N = length(P_C[:,1])                            # Number of generators
F = rand(100:500,N)                             # Random prod. prices
S = F*100                                       # Start up costs
T_up = rand(3:10,N)                             # Minimum operation time
T_down = rand(3:5,N)                            # Minimum shut down time
atol = 10e-6                                    # tolerance variable

# ========================= FUNCTIONS =========================
# Recursive functions for minimum operating and shut down time
function consecutive(T,n,t,sum=0)
    if T > 0
        sum += x[n,t + T]
        consecutive(T - 1,n,t,sum)
    else
        return sum
    end
end
function consecutive_down(T,n,t,sum=0)
    if T > 0
        sum += x[n,t + T] - 1
        consecutive(T - 1,n,t,sum)
    else
        return sum
    end
end

# ========================= MODEL =========================
m = Model(CPLEX.Optimizer)                                  # Model
@variable(m, x[1:N,1:T], Bin,start=0)                       # Unit activation
@variable(m, x_neg[1:N,1:T], Bin,start=0)                   # Negated unit activation
@variable(m, x_up[1:N,1:T], Bin)                            # Unit activation
@variable(m, x_down[1:N,1:T], Bin)                          # Unit activation
@variable(m, P_G[i=1:N,1:T])                                # Power generation

for t in 1:T                                                # Load balance
    @constraint(m, sum(P_G[:,t]) == P_D[t])
end


for n in 1:N                                                # Unit and time dependent constraints                                              
    for t in 1:T
        # Unit generation limit
        @constraint(m, P_C[n,1]*x[n,t] <= P_G[n,t])
        @constraint(m,P_G[n,t] <= P_C[n,2]*x[n,t])
        
        # start up and shut down, Indicator constraint for Bins
        if t > 1
            @constraint(m, !x_up[n,t]  => {x[n,t-1] - x[n,t] >= 0-atol}) # up is 0 when -1
            @constraint(m, !x_down[n,t] => {x[n,t-1] - x[n,t] <= 0+atol}) # down is 0 when 1
        else
            @constraint(m, !x_up[n,t]  => {0 - x[n,1] >= 0-atol})
        end
        
        # Minimum operation time, Indicator constraint for Bins
        if t > T - T_up[n]
            @constraint(m, x_up[n,t] => {consecutive(T - t,n,t) == T - t})
        else
            @constraint(m, x_up[n,t] => {consecutive(T_up[n],n,t) == T_up[n]})
        end
        
        # NOT IMPLEMENTED: Minimum cool down time, Indicator constraint for Bins
        if t > T - T_down[n]
            #@constraint(m, x_down[n,t] => {consecutive(T_down[n],n,t) == 0})
        else
            #@constraint(m, x_down[n,t] => {consecutive_down(T_down[n],n,t) == -T_down[n]})
        end
        #@constraint(m, !x_neg[n,t] => {x[n,t] <= 0+atol})
    end
end

@objective(m,Min,
    sum(P_G[:,1:T].*F[1:N])                         # Fuel cost
    + sum(x_up.*S[1:N])                             # Start-up cost
)

optimize!(m)                                                # Solve

# ========================= PLOT =========================
clrs = [round(Int,i) == 1 ? "green" : "red" for i in value.(x)]        # https://stackoverflow.com/questions/39790031/does-julia-have-a-ternary-conditional-operator
clrs1 = [round(Int,i) > 0 ? "grey" : "white" for i in value.(x_up) + value.(x_down)]        # https://stackoverflow.com/questions/39790031/does-julia-have-a-ternary-conditional-operator
p1 = plot(                                                      # Unit activation
    1:T,value.(P_G[:,1:T])',
    xtick=0:10:T,
    ylab = L"P_{unit} [MW]",ylim=(minimum(P_C),maximum(P_C)))
p2 = scatter(                                                   # Unit activation
    repeat(1:T, 1,N)',
    repeat(1:N, 1,N),
    color=clrs,
    ytick=(1:N, ["$(i), $(F[i]) \$, [$(P_C[i,1])-$(P_C[i,2])] MW" for i in 1:N]),
    legend=false,
    xtick=0:10:T,
    yflip=true,
    ylab = "unit")
scatter!(                                                   # Unit activation
    repeat(1:T, 1,N)',
    repeat(1:N, 1,N),
    color=clrs1,
    markershape = :cross,
    markerstrokealpha = 0.2,
    markerstrokewidth = 0.2,
    markersize = 3,
    ytick=(1:N, ["$(i), $(F[i]) \$, [$(P_C[i,1])-$(P_C[i,2])] MW" for i in 1:N]),
    legend=false,
    xtick=0:10:T,
    yflip=true)

p3 = plot(1:T,
    value.(P_G[:,1:T])'*F[1:N]  + (value.(x_up)[:,1:T]'*S[1:N]),
    ylab = "\$",
    xtick=0:10:T,
    xlab = L"t")
plt = plot(p1, p2, p3, layout = (3, 1), legend = false,size=(900,600))


@show plt
