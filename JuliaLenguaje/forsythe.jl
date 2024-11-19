using DifferentialEquations
using Plots

# Definimos los parámetros del problema
μ = 1 / 82.45
μ_star = 1 - μ

# Función que describe el sistema de ecuaciones diferenciales
function three_body!(du, u, p, t)
    x, y, vx, vy = u
    r1 = ((x + μ)^2 + y^2)^(3/2)
    r2 = ((x - μ_star)^2 + y^2)^(3/2)
    du[1] = vx
    du[2] = vy
    du[3] = 2 * vy + x - ((μ_star * (x + μ)) / r1) - ((μ * (x - μ_star)) / r2)
    du[4] = -2 * vx + y - ((μ_star * y) / r1) - ((μ * y) / r2)
end

# Condiciones iniciales y período de integración
u0 = [1.2, 0.0, 0.0, -1.04935751]
tspan = (0.0, 6.19216933)
prob = ODEProblem(three_body!, u0, tspan)

# Resolver el problema de ODE
sol = solve(prob, Tsit5(), abstol=1e-12, reltol=1e-5)

# Calcular la distancia mínima a la Tierra en la órbita
distances = [sqrt((u[1] + μ)^2 + u[2]^2) for u in sol.u]
min_distance = minimum(distances)

# Convertir la distancia mínima a millas (unidad Tierra-Luna = 238,000 millas)
min_distance_miles = min_distance * 238000

println("La distancia mínima del Apollo a la superficie de la Tierra es de: $(min_distance_miles - 4000) millas.")

# Graficar la órbita
plot(sol, vars=(1, 2), title="Órbita del tercer cuerpo", xlabel="x", ylabel="y", legend=false)
