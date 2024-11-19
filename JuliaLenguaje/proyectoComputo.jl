# Importamos el paquete DifferentialEquations.jl para trabajar con ecuaciones diferenciales
using DifferentialEquations

# Definimos los parámetros del problema, μ es la razón de las masas de la luna a la tierra
μ = 1 / 82.45  # Constante dada en el problema
μ_star = 1 - μ  # Complemento de μ, para simplificar las fórmulas

# Función que describe el sistema de ecuaciones diferenciales de primer orden
# du: Vector de derivadas
# u: Vector de variables de estado [x, y, vx, vy]
# p: Parámetros (no se usan en este caso)
# t: Tiempo
function three_body!(du, u, p, t)
    # Desempaquetamos las variables de estado
    x, y, vx, vy = u
    
    # Calculamos r1 y r2, las distancias desde el tercer cuerpo a los otros dos cuerpos
    r1 = ((x + μ)^2 + y^2)^(3/2)  # Distancia al cuerpo en (-μ, 0)
    r2 = ((x - μ_star)^2 + y^2)^(3/2)  # Distancia al cuerpo en (1 - μ, 0)
    
    # Ecuaciones diferenciales del problema
    du[1] = vx  # dx/dt = velocidad en x
    du[2] = vy  # dy/dt = velocidad en y
    du[3] = 2 * vy + x - ((μ_star * (x + μ)) / r1) - ((μ * (x - μ_star)) / r2)  # dvx/dt
    du[4] = -2 * vx + y - ((μ_star * y) / r1) - ((μ * y) / r2)  # dvy/dt
end

# Condiciones iniciales para las variables de estado [x(0), y(0), x'(0), y'(0)]
u0 = [1.2, 0.0, 0.0, -1.04935751]  # Dadas en el problema

# Intervalo de tiempo para la integración (inicio y fin)
tspan = (0.0, 6.19216933)  # Período del problema

# Definimos el problema de ODE
prob = ODEProblem(three_body!, u0, tspan)

# Resolvemos el problema usando el método de Runge-Kutta-Fehlberg de orden 4/5 (RK45)S
sol = solve(prob, Tsit5())

# Importamos Plots para visualizar la solución
using Plots

# Graficamos la solución en el plano xy (trayectoria del tercer cuerpo)
plot(sol, vars=(1, 2), title="Órbita del tercer cuerpo", xlabel="x", ylabel="y")
