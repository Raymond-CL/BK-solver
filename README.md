# BK-solver

A program written in Fortran that generates a 2-dimensional table of the dipole amplitude $N(r,Y)$ by solving the Balitsky-Kovchegov (BK) equation given an intial condition at $N(r,Y=0)$. The BK equation in coordinate (configuration) space is given as:

$$
\frac{\partial}{\partial Y} N(r_{01},Y) = \int d^2r_2 ~ {\rm ker}(r_{01},r_{02},r_{12}) \left[N(r_{02},Y)+N(r_{12},Y)-N(r_{01},Y)-N(r_{02},Y)N(r_{12},Y)\right]
$$

where the quark anti-quark ( $q\bar{q}$ ) dipole has transverse coordinate at $r_0$ and $r_1$ with size (separation distance) $r_{01}\equiv |\vec{r_0}-\vec{r_1}|$. The equation gives the evolution of this dipole amplitude $N(r,Y)$ along the rapidity $Y=\ln(s/s_0)=\ln(x_0/x)$ or equivalently along the momentum fraction $x$. The RHS of the equation performs an integration on $r_2$, which is the transverse coordinate of the radiated gluon ( $g$ ) by the dipole. In the large $N_c$ approximation, the gluon can be viewed as a $q\bar{q}$ pair, which gives another two dipoles of size $r_{02}$ and $r_{12}$. The integration is then simply the Riemann sum on all possible transverse coordinates of $r_2$ relative to $r_{01}$. ${\rm ker}$ is the kernel of the integration.

Most of the program functions can be configured by setting the `input.dat` file in root directory.

### input description:

- `rmin,rmax,rn` sets the range of the dipole dize in the table
  - 
