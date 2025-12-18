# From ChatGPT
import sympy as sp

# symbols
t = sp.symbols('t')
b0, b1, a0, a_star, tau_P = sp.symbols(
    'b0 b1 a0 a_star tau_P', positive=True
)

# function
a = sp.Function('a')(t)

# RHS
rhs = (
    (b0 + b1*(a0 - a))
    * sp.log((b0 + b1*(a0 - a)) / (b0 + b1*(a0 - a_star)))
    / (b1 * tau_P)
)

# ODE
ode = sp.Eq(sp.diff(a, t), rhs)

ode

sol = sp.dsolve(ode)

sol

# Constant of integration
a_init, C1 = sp.symbols('a_init C1')

a_expr = -a0*sp.exp(sp.exp(-C1*b1 - t/tau_P)) + a0 + a_star*sp.exp(sp.exp(-C1*b1 - t/tau_P)) - b0*sp.exp(sp.exp(-C1*b1 - t/tau_P))/b1 + b0/b1

C1_sol = sp.solve(sp.Eq(a_expr.subs(t, 0), a_init), C1)[0]
sp.simplify(C1_sol)

a_expr.subs(C1, C1_sol).simplify()

# (-b0*exp(exp(-t/tau_P)*log((a0*b1 - a_init*b1 + b0)/(a0*b1 - a_star*b1 + b0))) + b0 + b1*(-a0*exp(exp(-t/tau_P)*log((a0*b1 - a_init*b1 + b0)/(a0*b1 - a_star*b1 + b0))) + a0 + a_star*exp(exp(-t/tau_P)*log((a0*b1 - a_init*b1 + b0)/(a0*b1 - a_star*b1 + b0)))))/b1
