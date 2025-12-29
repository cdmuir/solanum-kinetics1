# Derivations for power function model

# Derive da / dP from a(P)

import sympy as sp

# symbols and unknown function
P = sp.symbols('P', real=True)
k, a_max, y0 = sp.symbols('k a_max y0', real=True)
a = sp.Function('a')

# ODE: da/dP = (k/(1+P)) * (a - a_max)
ode = sp.Eq(sp.diff(a(P), P), (k/(1 + P)) * (a(P) - a_max))

# solve with initial condition a(0) = y0
sol = sp.dsolve(ode, ics={a(0): y0})
print(sol)

# simplify to a power form: exp(k*log(1+P)) -> (1+P)**k
sol_simplified = sp.Eq(sol.lhs, sp.simplify(sol.rhs.rewrite(sp.Pow)))
print(sol_simplified)

# Rearrange a(P) to get P(a)

# symbols
P, a = sp.symbols('P a', real=True)
a_max, a0, k = sp.symbols('a_max a0 k', real=True)

# original equation
expr = sp.Eq(
    a,
    a_max + (P + 1)**k * (a0 - a_max)
)

# solve for P
sol = sp.solve(expr, P)
sol

# Get P* - P in terms of a* - a
# symbols
P_star, P, a_star, a = sp.symbols('P_star P a_star, a', real=True)

sp.simplify((((a_star - a_max)/(a0 - a_max))**(1/k) - 1) - (((a - a_max)/(a0 - a_max))**(1/k) - 1))

# Get k * (a - a_max) / (1 + P) in terms of a

sp.simplify(
  k * (a - a_max) / (1 + (((a - a_max)/(a0 - a_max))**(1/k) - 1))
)

# Simplify the whole da / dt expression
tau_P = sp.symbols('tau_P', real=True)
sp.simplify(((-(a_max - a_star)/(a0 - a_max))**(1/k) - ((a - a_max)/(a0 - a_max))**(1/k)) / tau_P * k*(a - a_max)/((a - a_max)/(a0 - a_max))**(1/k))

## Get a(t) straight from P(t)
import sympy as sp

# symbols
P, P_init, P_star, t = sp.symbols('P P_init P_star t', real=True)
a_init, a_max, a_star, a0, k, tau_P = sp.symbols('a_init a_max a_star a0 k tau_P', real=True)

# define P(a) and a(P)
P = lambda a: ((a - a_max)/(a0 - a_max))**(1/k) - 1
a = lambda P: a_max + (P + 1)**k*(-a_max + a0)

P_init = P(a_init)
P_star = P(a_star)

P_t = P_star - (P_star - P_init) * sp.exp(-t / tau_P)

P_t_simplified = sp.simplify(P_t)
a(P_t_simplified).simplify()
