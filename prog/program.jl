using Plots
using Printf
using Polynomials
using QuadGK

# WSZYSTKIE FUNKCJE, KTÓRYCH UŻYŁEM DO WYKONANIA OBLICZEŃ


# funkcja wyznaczająca równoodległe węzły z przedziału [a, b]

function nodes(a, b, n) 
    t = zeros(n + 1)

    h = (b - a) / n

    for i = 0:n
        t[i + 1] = a + h * i
    end
    
    return t
end


# funkcja wyznaczająca wartości drugich pochodnych funkcji sklejanej 3 stopnia

function secDerValues(f, a, b, n, t)
    h = (b - a) /  n
    
    u = zeros(n-1)
    b = zeros(n)
    v = zeros(n-1)
    lambdas = zeros(n + 1)
    lambdas[1] = 0
    lambdas[n+1] = 0

    for i=1:n
        b[i] = 6(f(t[i + 1]) - f(t[i])) / h
    end
    
    u[1] = 2 * (h + h) 
    v[1] = b[2] - b[1]
    
    for i = 2:n-1
        u[i] = 2 * (h + h) - h * h / u[i - 1]
        v[i] = b[i+1] - b[i] - h * v[i - 1] / u[i - 1]
    end
    
    for i = n:-1:2
        lambdas[i] = (v[i - 1] - h * lambdas[i + 1]) / u[i - 1]
    end

    return lambdas
end


# funkcja wyznaczająca współczynniki każdego wielomianu 3 stopnia funkcji sklejanej

function getSplineGraph(f, a, b, n) 
    h = (b - a) / n
    t = nodes(a, b, n)
    l = secDerValues(f, a, b, n, t)
    y = zeros(n + 1)
    coeff = zeros(n, 4)
    
    for i=1:n+1
        y[i] = f(t[i])
    end
    
    for i=1:n
        coeff[i, 1] = (-l[i] / (6*h) + l[i+1] / (6*h))
        coeff[i, 2] = (l[i] / (2*h) * t[i + 1] - l[i + 1] / (2 * h) * t[i])
        coeff[i, 3] = (-l[i] / (2 * h) * (t[i + 1])^2 + l[i + 1] / (2 * h) * (t[i])^2 + y[i + 1] / h - 
            l[i + 1] * h / 6 - y[i] / h + l[i] * h / 6)
        coeff[i, 4] = (l[i] / (6 * h) * (t[i+1])^3 - l[i + 1] / (6 * h) * (t[i])^3 - y[i+1] / h * t[i] + 
            l[i+1] * h / 6 * t[i] + y[i] / h * t[i + 1] - l[i] * h / 6 * t[i + 1])
    end

    return (coeff)
end


# funkcja zwracająca spróbkowane wartości funkcji sklejanej w punktach xs

function getSplineValues(a, b, n, coeff, numberOfPoints = 100) 
    
    xs = LinRange(a, b, numberOfPoints * n)
    ys = zeros(numberOfPoints * n)
    
    for i=1:n 
        g(x) = coeff[i, 1] * x^3 + coeff[i, 2] * x^2 + coeff[i, 3] * x + coeff[i, 4]

        for j=1:numberOfPoints
            ys[(i - 1) * numberOfPoints + j] = g(xs[(i - 1) * numberOfPoints + j])
        end 
    end
    
    return xs, ys
end


# funkcja zwracająca wartośc całki funkcji sklejanej dla zadanych wspólczynników

function integrateSpline(f, a, b, n, x)
    h = (b - a) / n
    t = nodes(a, b, n)
    
    coeff = getSplineGraph(f, a, b, n)
    
    lastIndex = 1
    
    while x >= t[lastIndex]
        lastIndex = lastIndex + 1
    end
    
    lastIndex = lastIndex - 1
    sum = 0
    
    for i=1:lastIndex
        g(x) = coeff[i, 1] * x^4 / 4 + coeff[i, 2] * x^3 / 3 + coeff[i, 3] * x^2 / 2 + coeff[i, 4] * x
        area = g(t[i+1]) - g(t[i])
        
        sum = sum + area
    end
    
    return sum
end


# funkcja zwracająca błąd bezwzględny z wartości całki funkcji sklejanej

function compare(f, a, b, n, x)
    int = quadgk(f, a, x)
    intSpline = integrateSpline(f, a, b, n, x)
    
    @printf("Błąd bezwzględny dla x = %f i n = %f:\n", x, n)
    @printf("%d & %f & %f \\\\", n, intSpline, abs(int[1] - intSpline))
    return abs(int[1] - intSpline), intSpline
end


# funkcja rysująca oba wykresy, funkcji przybliżanej i funkcji sklejanej

function drawSpline(f, a, b, n)
    coeff = getSplineGraph(f, a, b, n)
    g = getSplineValues(a, b, n, coeff)

    default(tickfont = (8, :red), framestyle = :zerolines)

    plot(g[1], g[2], linewidth = 2, label = "Funkcja sklejana", legend = :bottomright)
    plot!(f, a, b, linewidth = 2, label = "Funkcja przybliżana") 
end