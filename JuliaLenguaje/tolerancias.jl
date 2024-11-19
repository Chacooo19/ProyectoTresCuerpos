using DifferentialEquations

# Definimos los parámetros del problema
μ = 1 / 82.45
μ_star = 1 - μ

# Función que describe el sistema de ecuaciones diferenciales
function three_body!(du, u, p, t)
    x, y, vx, vy = u
    r1 = ((x + μ)^2 + y^2)^(3/2)
    r2 = ((x - μ_star)^2 + y^2)^(3/2)
    
    du[1] = vx  # dx/dt = vx
    du[2] = vy  # dy/dt = vy
    du[3] = 2 * vy + x - ((μ_star * (x + μ)) / r1) - ((μ * (x - μ_star)) / r2)
    du[4] = -2 * vx + y - ((μ_star * y) / r1) - ((μ * y) / r2)
end

# Condiciones iniciales
u0 = [1.2, 0.0, 0.0, -1.04935750983031990726]
tspan = (0.0, 6.19216933131963970674)

# Definir y resolver el problema con tolerancias específicas
prob = ODEProblem(three_body!, u0, tspan)
sol = solve(prob, Tsit5(), abstol=1e-9, reltol=1e-6)  # Tolerancias recomendadas

# Graficar la solución
using Plots
plot(sol, vars=(1, 2), title="Órbita del tercer cuerpo con tolerancias", xlabel="x", ylabel="y")
