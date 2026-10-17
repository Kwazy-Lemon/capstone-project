# Linear vs. Nonlinear Components in MNA

Notes on how circuit components are handled in the simulator, and why nonlinear devices (starting with the diode) require Newton's method.

## Contents

1. [Linear components](#1-linear-components)
2. [Nonlinear components: the diode](#2-nonlinear-components-the-diode)
3. [Newton's method](#3-newtons-method)
4. [Applying Newton's method to a diode](#4-applying-newtons-method-to-a-diode)
5. [Companion model](#5-companion-model)
6. [Implementation outline](#6-implementation-outline)

---

## 1. Linear components

### Resistor

$$i = G\,v$$

where $G = 1/R$ is the conductance. Stamping every component into the MNA matrix gives a single linear system:

$$G\,x = b$$

The matrix entries do not depend on the solution, so one linear solve is enough.

### Capacitor and inductor

$$i_C = C\,\frac{dv}{dt}$$

Derivative terms are handled by **Backward Euler**, which turns each time step into another linear solve. For the MNA system $G\,x(t) + C\,\dot{x}(t) = b(t)$, approximate

$$\dot{x} \approx \frac{x_{n+1} - x_n}{h}$$

where $h$ is the time step and $x_n$ is already known. Substituting gives

$$\left(G + \frac{C}{h}\right) x_{n+1} = b_{n+1} + \frac{C}{h}\,x_n$$

This is a plain $A\,x = b$ problem at every time step.

---

## 2. Nonlinear components: the diode

A diode's current follows the Shockley equation:

$$i_d = I_s\left(e^{v_d/(n V_t)} - 1\right)$$

Two issues arise from this equation:

1. The slope of $i$ vs. $v$ is not constant. For example, at 0.3 V it is tiny and at 0.7 V it is enormous.
2. To stamp the diode into the matrix $G$ its conductance is needed, but that conductance depends on $v_d$, which is part of the answer that has not been found yet.

In other words, the solution is needed to build the matrix, and the matrix is needed to get the solution. The system is therefore no longer

$$G\,x = b$$

but becomes

$$F(x) = 0$$

where $F$ is a nonlinear function of the unknown vector $x$. It cannot be solved with a single linear solve.

For the linear case, $F$ is simply the residual of the original system:

$$F(x) = G\,x - b$$

In general:

$$F(x) = G_{lin}\,x + f_{nl}(x) - b = 0$$

where $f_{nl}(x)$ contains the nonlinear device currents.

---

## 3. Newton's method

Newton's method solves nonlinear equations by repeatedly replacing them with a linear approximation.

### One variable

To find $x$ such that $f(x) = 0$:

1. Guess $x_0$.
2. Draw the tangent line at $x_0$. The tangent has slope $f'(x_0)$.
3. Find where the tangent crosses zero using point-slope form:

$$x_1 = x_0 - \frac{f(x_0)}{f'(x_0)}$$

This is the next guess.

4. Repeat until the change is tiny.

### Multiple unknowns

With multiple unknowns the derivative becomes the **Jacobian matrix** $J$:

$$J_{i,j} = \frac{\partial F_i}{\partial x_j}$$

$J_{i,j}$ answers: "if node $j$'s voltage is nudged slightly, how much does the imbalance in equation $i$ change?"

Near the current guess, $F(x + \Delta x) \approx F(x) + J\,\Delta x$. Setting this to zero gives a linear system for the correction $\Delta x$:

$$J(x_k)\,\Delta x = -F(x_k), \qquad x_{k+1} = x_k + \Delta x$$

### Each Newton iteration

1. Evaluate $F$
2. Evaluate $J$
3. Do a linear solve for $\Delta x$
4. Update $x$ and check convergence

---

## 4. Applying Newton's method to a diode

With Modified Nodal Analysis, the circuit equations are KCL at each node: the sum of currents leaving the node equals 0. $F(x)$ is the KCL residual.

- The **linear parts** contribute $G_{lin}\,x - b$ to $F$, and $G_{lin}$ to $J$.
- The **diode** (anode node $a$, cathode node $k$, $v_d = v_a - v_k$) contributes its current

$$i_d = I_s\left(e^{v_d/(n V_t)} - 1\right)$$

to $F$. This current leaves node $a$ and enters node $k$.

- Its derivative, the small-signal conductance, goes into $J$:

$$g_d = \frac{di_d}{dv_d} = \frac{I_s}{n V_t}\,e^{v_d/(n V_t)}$$

The conductance is stamped using the same four-entry pattern as a resistor:

|       | $a$    | $k$    |
|-------|--------|--------|
| **$a$** | $+g_d$ | $-g_d$ |
| **$k$** | $-g_d$ | $+g_d$ |

At every iteration the diode behaves like a resistor of value $1/g_d$, and that value changes each iteration as $v_d$ changes.

### Algorithm

1. Guess $x_0$ (default: 0)
2. Compute $F(x_k)$ and $J(x_k)$
3. Solve $J\,\Delta x = -F(x_k)$
4. Update $x_{k+1} = x_k + \Delta x$
5. Repeat until the change is tiny

### Equivalent form: solving directly for $x_{k+1}$

Substituting $\Delta x = x_{k+1} - x_k$:

$$J\,x_{k+1} = J\,x_k - F(x_k)$$

This solves directly for the new $x_{k+1}$ instead of the correction. The right side is known and the left side is a matrix times the unknown, so it is a general linear system $A\,x = b$.

---

## 5. Companion model

The companion model is the same Newton iteration, rearranged so it reuses the existing linear stamping code.

Linearize the diode current around the current guess $v_k$:

$$i(v) \approx i_d(v_k) + g_d\,(v - v_k) = g_d\,v + \underbrace{\left[i_d(v_k) - g_d\,v_k\right]}_{I_{eq}}$$

Reading the result:

- $g_d\,v$ is a **resistor** (conductance $g_d$)
- $I_{eq} = i_d(v_k) - g_d\,v_k$ is a **constant current source**

A conductance $g_d$ in parallel with a current source $I_{eq}$ is a circuit that can be stamped with the existing linear component code.

### At each iteration

1. Compute $g_d$ and $I_{eq}$ from the current guess of the diode voltage
2. Stamp them like a resistor and a current source
3. Solve the ordinary linear system for the new voltages
4. Check whether the voltages changed by less than a tolerance. If not, repeat from step 1

### Newton form vs. companion form

| Newton form | Companion form |
|---|---|
| Build $J$ and $F$, solve for $\Delta x$, then add | Stamp $g_d$ and $I_{eq}$, solve for $x_{new}$ directly |
| Needs a residual function and a Jacobian function | Needs one "stamp at this guess" function |

Both produce the same iterates.

---

## 6. Implementation outline

```text
x = initial guess (zeros)
loop:
    build matrix from linear stamps          # constant, can be built once
    add nonlinear stamps evaluated at x      # rebuilt every iteration
    solve  ->  x_new
    if change is small (and residual is small):
        done
    x = x_new
```

Practical issues to handle:

- **Exponential overflow:** a large trial $v_d$ makes $e^{v_d/(nV_t)}$ overflow. Limit the voltage step per iteration (SPICE uses `pnjlim`).
- **Convergence is not guaranteed:** poor initial guesses can oscillate or diverge. Fallbacks include source stepping and gmin stepping.
- **Stopping criterion:** check both $\lVert \Delta x \rVert$ and $\lVert F \rVert$ against tolerances (abstol/reltol).
- **Assembler structure:** keep constant (linear) stamps separate from per-iteration (nonlinear) stamps.
- **Transient analysis:** a Newton loop runs inside every time step.
- **Component interface:** each nonlinear component provides a method that, given the current guess, returns either its $(F, J)$ contribution or the equivalent $(g, I_{eq})$ stamp.
