# BK-solver

>_*_This program is currently under construction!!!_*_

A program written in Fortran that generates a 2-dimensional table of the dipole amplitude $N(r,Y)$ by solving the Balitsky-Kovchegov (BK) equation given an intial condition at $N(r,Y=0)$. The BK equation in coordinate (configuration) space is given as:

$$
\frac{\partial}{\partial Y} N(r_{01},Y) = \int d^2\vec{r_2} ~ {\rm ker}(r_{01},r_{02},r_{12}) \left[N(r_{02},Y)+N(r_{12},Y)-N(r_{01},Y)-N(r_{02},Y)N(r_{12},Y)\right]
$$

where the quark anti-quark ( $q\bar{q}$ ) dipole has transverse coordinate at $r_0$ and $r_1$ with size (separation distance) $r_{01}\equiv |\vec{r_0}-\vec{r_1}|$. The equation gives the evolution of this dipole amplitude $N(r,Y)$ along the rapidity $Y=\ln(s/s_0)=\ln(x_0/x)$ or equivalently along the momentum fraction $x$. The RHS of the equation performs an integration on $r_2$, which is the transverse coordinate of the radiated gluon ( $g$ ) by the dipole. In the large $N_c$ approximation, the gluon can be viewed as a $q\bar{q}$ pair, which gives another two dipoles of size $r_{02}$ and $r_{12}$. The integration is then simply the Riemann sum on all possible transverse coordinates of $r_2$ relative to $r_{01}$. ${\rm ker}$ is the kernel of the integration.

## what the program *is* and *is not* for

- The dipole amplitude $N(r,Y)$ in configuration space is usually used to calculate cross-sections of particle productions in high energy nuclear collisions especially in the saturation (small- $x$ ) regime. It needs to be Fourier transformed into momentum space:

$$
\tilde{N}(k,x) = \int d^2\vec{r} ~ e^{-i\vec{r}\cdot\vec{k}} \left[1-N(|\vec{r}|,Y=\ln(x_0/x))\right]
$$

where $\tilde{N}$ is interpret as the unintegrated gluon distribution.
- The program calculates the dipole amplitude in coordinate space ( $\vec{r}$ ), not in momentum space ( $\vec{k}$ ). I do plan on writing a program directly in momentum space, but kinematic constrain will be difficult to apply.
- The program also does not include impact parameter ( $\vec{b}$ ) dependence. This is actually a very important part that contributes in heavy-ion collisions, but the overall procedure is complicated and will be included in the next iteration of the program.

## getting started

Most of the programs that I wrote is completely *stand-alone*, but if you are a fresh beginner who just got into this field, and just so happens to find this program in your github search result. Well, you can follow the procedures below to get the program to run.

- You need a linux environment and have git. Clone this repository to your local directory.
- You'll need fortran compiler in your system. On WSL or Ubuntu:
```
sudo apt install gfortran build-essential makedepf90
```
which installs the necessary compilers.
- in the root directory, the command `make` should compile all the codes and produce an executable. While `make run` should produce and run the executable. And `make clean` deletes all unnecessary by-products of the program without deleting the source codes.

## input description:

Most of the program functions can be configured by setting the [input.dat](./input.dat) file in root directory.
Some of the more subtle 'user defined' modifications can be set in the source code [BK.f90](./src/BK.f90) file, but the program needs to be recompile for the modifications to take place.

- `rmin,rmax,rn` sets the range of the dipole dize in the table
  - theoretically the dipole size ranges from 0 to infinity. But because $r\rightarrow0$ is singular and the table arrange the value in logarithmic scale, `rmin` should be stricly greater than 0.
  - the range of $10^{-10} < r < 10^{+5}$ should cover most phase space. If you want to evolve to a larger $Y$ value, then you should extend `rmin` to smaller value.
- `ymin,ymax,yn` set the range of rapidity (range of evolution)
  - since the initial condition is set to $Y=0$ (mostly), `ymin` should be set to 0.
  - `ymax` at 60 corresponds to momentum fraction $x\sim 10^{-30}$, enough for current particle colliders.
- `IniCnd` sets the initial condition of the dipole amplitude at $N(r,Y=0)$
  - `1` GBW- $\gamma$ prescription
  - `2` MV- $\gamma$ prescription
  - `3` user defined
- `EvoMth` set the evolution method (method of solving differential equation)
  - Given $N(r,Y)$, evolve to next step $N(r,Y+h)$ with step width $h$ using the standard Runge-Kutta method at different order
  - `1` 1st-order Runge-Kutta at $\mathcal{O}(h^1)$ (Euler's method)
  - `2` 2nd-order Runge-Kutta at $\mathcal{O}(h^2)$
  - `3` 4th-order Runge-Kutta at $\mathcal{O}(h^4)$ (RK4)
- `IntMth` sets the integration method of $d^2r_2$
  - `1` $d^2r_2=dx~dy$, generates cartesian coordinates $(x,y)$ that is used to determine $r_{02}$ and $r_{12}$, with symmetric factor = 4.
  - `2` $d^2r_2=d\phi~rdr$, generates polar coordinates $(r,\phi)$ relative to the dipole center and determine $r_{02}$ and $r_{12}$, with symmetric factor = 4.
  - `3` $d^2r_2=d\theta~r_{02}dr_{02}$, generates polar coordinates $(r_{02},\theta)$ relative to $r_0$ and determine $r_{12}$, with symmetric factor = 2.
- `EvoKer` sets the BK evolution kernel ${\rm ker}$
  - `1` leading-order (parent dipole) prescription
  - `2` Balitsky prescription
  - `3` Kovchegov-Weigert prescription
- `RunCup` sets the running coupling in the kernel
  - `1` fixed coupling (default is 0.2, can modify in `BK.f90`)
  - `2` running with $r$
  - `3` user defined
- `IntPol` sets the interpolation method when reading the $N(r,Y)$ table for integration
  - `1` $N(r)$ vs $r$ is assumed to be linear
  - `2` $\ln(N(r))$ vs $\ln(r)$ is assumed to be linear
  - `3` natural cubic spline in the `./nr/` package (super slow)
- `ncall,itmax` sets the vegas integration input.
  - *vegas* is an adaptive multi-dimensional Monte-Carlo integration routine.
  - `ncall` sets the number of samples per iteration
  - `itmax` sets the maximum iteration number
  - With every iteration, vegas adapts to the behaviour of the integrand and sets the sample distribution accordingly, this is the *warm-up stage* of vegas. Then the *final stage* is performed by using minimum number of iteration and maximum number of samples to increase the precision of the integration.
  - In the *warm-up stage*, 10000 samples by `itmax` iterations are run. In the *final stage*, `ncall` samples by 1 iterations are run.
  
## program improvements and side notes

- One of the main improvements that could be implemented in this program is the inclusion of the impact parameter $b$ dependence. However, this feature is very complicated and would result in a separate version of the program.
- Parallelization of the codes, in Fortran. Although the new Fortran standard has openmp features, it is still difficult to implement it in the NR package. But it can be done.
- Momemtum space? need to solve the kinematic constrain problem first.
- NLO? logarithmic divergences expected. Resummation needed.
- let user choose whether to include the non-linear term or not. The $-N(r_{02})N(r_{12})$ term that prevents $N(r)$ from exceeding 1 (unitary bound). (BFKL or BK)
- For sufficiently large $r$, $N(r,Y)\rightarrow1$, then $\partial/\partial Y~N(r,Y)\rightarrow 0$ (BK).
- $r_2$ cannot be close to $r_0$ or $r_1$, which makes either $r_{02}$ or $r_{12}$ very small, resulting in divergence in the kernel.
- this divergence is then cancelled be the smallness of $N(r)$
