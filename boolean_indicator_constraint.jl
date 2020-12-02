using JuMP, Cbc             

# ========================= DATA =========================
P_C  = [50 200; 25 250; 75 300; 100 400; 125 500; 150 600; 175 700; 200 800; 225 900; 100 1000;]
P_D = [3268.78 3158.69 3036.5 2946.15 2891.74 2909.41 3014.69 3104.31 3243.79 3547.45 3808.79 4134.67 4157.68 4048.36 4184.9 4052.32 4229.83 4510.66 4403.63 4158.34 3916.65 3669.91 3421.24 3164.72]
T = length(P_D)                                             # Number of time steps
N = length(P_C[:,1])                                        # Number of generators
F = rand(100:500,N)                                         # Random prod. prices
S = F*100
T_up = repeat([2],10,1)
S_min = repeat([3],10,1)

# ========================= MODEL =========================
next(t, i) = x[((t - 1 + i) % T) + 1]
m = Model(Cbc.Optimizer)                                    # Model
@variable(m, x[1:N,1:T], Bin,start=0)                       # Unit activation
@variable(m, P_G[i=1:N,1:T])                                # Power generation
@variable(m, U_up[i=1:N,1:T],)                              # Power generation
xx = @expression(m,[x[:,1] (x[:,2:T] - x[:,1:T-1])])
for t in 1:T                                                # Load balance
    @constraint(m, sum(P_G[:,t]) == P_D[t])
end
for n in 1:N                                                # Unit generation limit
    for t in 1:T
        @constraint(m, P_C[n,1]*x[n,t] <= P_G[n,t])
        @constraint(m,P_G[n,t] <= P_C[n,2]*x[n,t])
    end
end

@constraint(m, xxxx[n=1:N,t=2:T], (!x[n,t-1] && x[n,t]) => {next(t, 1) + next(t, 2) == 2})

@objective(m,Min,
    sum(P_G[:,1:T].*F[1:N])
)

optimize!(m)                                                # Solve
