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

# Condiciones iniciales y período de integración (máxima precisión de Shampine)
u0 = [1.2, 0.0, 0.0, -1.04935750983031990726]  # Incluyendo todos los decimales
tspan = (0.0, 6.19216933131963970674)  # Máxima precisión propuesta por Shampine

# Definimos el problema de ODE
prob = ODEProblem(three_body!, u0, tspan)

# Resolver el problema de ODE con alta precisión
sol = solve(prob, Vern9(), saveat=0.0001, abstol=1e-16, reltol=1e-14)

# Calcular la distancia mínima con mayor precisión
times = 0:0.0001:tspan[2]  # Rango de tiempo más fino
min_distance = minimum([sqrt((sol(t)[1] + μ)^2 + sol(t)[2]^2) for t in times])
min_distance_miles = min_distance * 238000 - 4000

println("La distancia mínima del Apollo a la superficie de la Tierra es de: $min_distance_miles millas.")

# Graficar la órbita
plot(sol, vars=(1, 2), title="Órbita del tercer cuerpo", xlabel="x", ylabel="y", legend=false)
