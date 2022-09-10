# BK-solver

>_*_This program is currently under construction!!!_*_

A program written in Fortran that generates a 2-dimensional table of the dipole amplitude $N(r,Y)$ by solving the Balitsky-Kovchegov (BK) equation given an intial condition at $N(r,Y=0)$. The BK equation in coordinate (configuration) space is given as:

$$
\frac{\partial}{\partial Y} N(r_{01},Y) = \int d^2\vec{r_2} ~ {\rm ker}(r_{01},r_{02},r_{12}) \left[N(r_{02},Y)+N(r_{12},Y)-N(r_{01},Y)-N(r_{02},Y)N(r_{12},Y)\right]
$$

where the quark anti-quark ( $q\bar{q}$ ) dipole has transverse coordinate at $r_0$ and $r_1$ with size (separation distance) $r_{01}\equiv |\vec{r_0}-\vec{r_1}|$. The equation gives the evolution of this dipole amplitude $N(r,Y)$ along the rapidity $Y=\ln(s/s_0)=\ln(x_0/x)$ or equivalently along the momentum fraction $x$. The RHS of the equation performs an integration on $r_2$, which is the transverse coordinate of the radiated gluon ( $g$ ) by the dipole. In the large $N_c$ approximation, the gluon can be viewed as a $q\bar{q}$ pair, which gives another two dipoles of size $r_{02}$ and $r_{12}$. The integration is then simply the Riemann sum on all possible transverse coordinates of $r_2$ relative to $r_{01}$. ${\rm ker}$ is the kernel of the integration.

Most of the program functions can be configured by setting the [input.dat](./input.dat) file in root directory.

## what the program *is* and *is not* for

- The dipole amplitude $N(r,Y)$ in configuration space is usually used to calculate cross-sections of particle productions in high energy nuclear collisions especially in the saturation (small- $x$ ) regime. It needs to be Fourier transformed into momentum space:

$$
\tilde{N}(k,x) = \int d^2\vec{r} ~ e^{-i\vec{r}\cdot\vec{k}} \left[1-N(|\vec{r}|,Y=\ln(x_0/x))\right]
$$

where $\tilde{N}$ is interpret as the unintegrated gluon distribution.
- The program calculates the dipole amplitude in coordinate space ( $\vec{r}$ ), not in momentum space ( $\vec{k}$ ). I do plan on writing a program directly in momentum space, but kinematic constrain will be difficult to apply.
- The program also does not include impact parameter ( $\vec{b}$ ) dependence. This is actually a very important part that contributes in heavy-ion collisions, but te overall procedure is complicated and will be included in the next iteration of the program.

## input description:

- `rmin,rmax,rn` sets the range of the dipole dize in the table
  - theoretically the dipole size ranges from 0 to infinity. But because $r\rightarrow0$ is singular and the table arrange the value in logarithmic scale, `rmin` should be stricly greater than 0.
  - the range of $10^{-10} < r < 10^{+5}$ should cover most phase space. If you want to evolve to a larger $Y$ value, then you should extend `rmin` to smaller value.
- `ymin,ymax,yn` set the range of rapidity (range of evolution)
  - since the initial condition is set to $Y=0$ (mostly), `ymin` should be set to 0.
  - `ymax` at 60 corresponds to momentum fraction $x\sim 10^{-30}$, enough for current particle colliders.
- `IniCnd` sets the initial condition of the dipole amplitude at $N(r,Y=0)$
  - 1. GBW- $\gamma$ prescription
  - 2. MV- $\gamma$ prescription
  - 3. user defined
- `EvoMth` set the evolution method (method of solving differential equation)
  - Given $N(r,Y)$, evolve to next step $N(r,Y+h)$ with step width $h$ using the standard Runge-Kutta method at different order
  - 1st-order Runge-Kutta at $\mathcal{O}(h^1)$ (Euler's method)
  - 2nd-order Runge-Kutta at $\mathcal{O}(h^2)$
  - 4th-order Runge-Kutta at $\mathcal{O}(h^4)$ (RK4)
- `IntMth` sets the integration method
