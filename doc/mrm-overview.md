# mRM - the Multiscale Routing Model <!--! {#mRM} -->

<!--! [TOC] -->

This section provides a short description of the method used in
mRM. For further details and an analysis over the European domain,
please refer to Thober et al. 2019.

<a name="findif"></a>

## Finite Difference Approximation of Kinematic Wave Equation

The multiscale Routing Model mRM uses the kinematic wave equation to
describe the water flow within a stream as

<a name="mrm_kinematic_wave"></a>

$$
\begin{equation}
\frac{\partial Q}{\partial t} + c \frac{\partial Q}{\partial x} = 0
\end{equation}
$$

where $Q$ ($\mathrm{m}^3\,\mathrm{s}^{-1}$) is streamflow, $x$ (m) the space dimension,
$t$ (s) the time dimension, and $c$ ($\mathrm{m}\,\mathrm{s}^{-1}$) the celerity. The
kinematic wave equation is a simplification of the Saint-Venant
equations. The derivation is based on the assumption
that the continuity equation is sufficient to describe the movement of
flood waves. In detail, a constant river bed slope and time-invariant
celerity $c$ have to be assumed. mRM employs
the classical finite difference weighted approximation on a four points
scheme to solve the equation above. A short derivation is as follows:

The partial derivatives within the above equation are represented as
finite differences between four values, that means on two locations at
two points in time:

$$
\begin{align*}
   \frac{\partial Q}{\partial t} & \approx \frac{\epsilon \left(Q(x_j,t_{i+1})
     - Q(x_j, t_i)\right) + (1 -
   \epsilon)\left(Q(x_{j+1}, t_{i+1})-Q(x_{j+1}, t_i)\right)}{\Delta t} , \\
   \frac{\partial Q}{\partial x} & \approx \frac{\varphi
     \left(Q(x_{j+1}, t_{i+1})
     - Q(x_j, t_{i+1})\right) + (1 -
   \varphi)\left(Q(x_{j+1}, t_i)-Q(x_j, t_i)\right)}{\Delta x} ,
\end{align*}
$$

where $j$ denotes the spatial location (i.e., reach id) and $i$ denotes
the timestep. $\epsilon$ is a space-weighting factor and $\varphi$ is a
time-weighting factor. mRM uses a rectangular grid to represent the
river network with a river reach in the model connecting two center grid
locations. Different spatial locations are separated by $\Delta x$ and
time steps by $\Delta t$. The two weighting factors, $\epsilon$ and
$\varphi$, can be chosen between 0 and 1, but the numerical solution
becomes unstable for $\epsilon >$ 0.5. The
numerical diffusion depends linearly on $\epsilon$,
with 0 implying full numerical diffusion and 0.5 no numerical diffusion,
respectively.

Setting $\varphi$ to 0.5, which represents a time-centered scheme, and
substituting the above equation into the
[kinematic wave equation](#mrm_kinematic_wave) results in the classical linear equation:

<a name="mrm_recurrence"></a>

$$
\begin{equation}
 Q(x_{j+1}, t_{i+1}) = C_1 Q(x_j, t_{i+1}) + C_2 Q(x_j, t_i) + C_3
 Q(x_{j+1}, t_i) ,
\end{equation}
$$

with the coefficients $C_1$, $C_2$ and $C_3$ being:
<a name="mrm_coefficients"></a>

$$
\begin{align}
   C_1 & = \frac{-2 \Delta x \epsilon + c \Delta t}{2 \Delta x (1 -
     \epsilon) + c \Delta t} , \nonumber \\
   C_2 & = \frac{2 \Delta x \epsilon + c \Delta t}{2 \Delta x (1 -
     \epsilon) + c \Delta t} , \\
   C_3 & = \frac{2 \Delta x (1 - \epsilon) - c \Delta t}{2 \Delta x (1 -
     \epsilon) + c \Delta t} . \nonumber
\end{align}
$$

The coefficients $C_1$, $C_2$ and $C_3$ add up to one. The spatial
resolution at which mRM is applied is called
"routing" resolution in the following.

<a name="slope_to_c"></a>

## Stream Celerity Parametrization based on Terrain Slope

Two parametrizations of the [coefficient expressions](#mrm_coefficients) are available in
mRM: firstly the regionalized Muskingum-Cunge (rMC) parametrization with
a fixed time step, and
secondly a parametrization using spatially varying
celerities in combination with an adaptive time step (aTS). A short
summary of the former is presented in Subsection [Regionalized Muskingum-Cunge parametrization](#rMC) and is referred to
as rMC in the following. The latter is described in detail in this and
the next section and is referred to as aTS scheme.

The aTS calculates stream celerity as a function of terrain slope using
a simple relationship:

$$
\begin{equation}
 c_i = \gamma \sqrt{s_i},
\end{equation}
$$

where $c_i$, $s_i$ and $\gamma$ are celerity, terrain slope and a free
model parameter, respectively, and $i$ is the grid cell index.

The above equation is applied at the resolution of the Digital
Elevation Model (DEM), from which terrain slope is derived. A median
absolute deviation (MAD) filter is applied to the
high resolution slope data along the path of the main river with a
threshold value of 2.25 to correct for outliers. Large slope outliers
happen easily in DEMs, for example, when the river flows in a valley
and one grid cell represents the valley while the next grid cell
represents the hill top. A minimum river bed slope of 0.1% is
further assumed. The celerities are then upscaled to the routing
resolution, i.e., the resolution at which the kinematic wave equation
is solved. The upscaling is by averaging with
the harmonic mean, the correct averaging operator for celerities (or
speed). The upscaling considers also only those high-resolution grid
cells that align with the path of the main river because aTS only
considers the flow in the main river reach, assuming that travel times
in the main reach dominate the routing process in
tributaries.

<a name="ada_time"></a>

## Adaptive Time Step (aTS) Implementation

The aTS scheme uses an adaptive time step to calculate the linear
coefficients in the [routing recurrence](#mrm_recurrence). The basic idea is that the
time step should be such that water has not been transported further
than $\Delta x$ during a single step. This condition is generally known
as Courant-Friedrich-Lewy criterium, which is a necessary condition for
numerical stability of finite difference schemes
\cite Courant:1928ch. The condition can be expressed as:

$$
\begin{equation}
 C_r = \frac{c \Delta t}{\Delta x} \le C_{max} = 1
\end{equation}
$$

where $C_{max}$ is the Courant number. aTS uses a Courant number
of 1. The Courant condition couples the applied spatial resolution
with the integration time step of the finite difference
scheme. Celerities $c_i$ are typically in the order of a few
$\mathrm{m}\,\mathrm{s}^{-1}$ The time step $\Delta t$ is chosen such
that it does not deviate too much from the Courant number
$C_{max}$ to keep computational
demand to a minimum. This implies that $\Delta t$ ranges from a
few minutes for high-resolution grids to a few hours for continental
scale applications.

In detail, aTS chooses $\Delta t$ from a set of prescribed values such
that $c_i\Delta t/\Delta x$ is close to but less than 1 for all
celerities $c_i$. The prescribed values range from one minute to one
day, namely: 1 min, 2 min, 3 min, 4 min, 5 min, 6 min, 10 min,
12 min, 15 min, 20 min, 30 min, 1 h, 2 h, 3 h, 4 h, 6 h, 8 h,
12 h, and 1 day. The choice of these values is motivated from the fact
that these represent multiples or dividers of hourly and daily time
steps. These time steps allow in principle model applications from
100 m to 100 km, for typical celerities around 1.5 $\mathrm{m}\,\mathrm{s}^{-1}$.

Note that the chosen time step depends only on the spatial resolution
and is independent of the time resolution of the provided forcing. For
example, applying aTS at 12 km resolution using a celerity of $c$~=
1.5 $\mathrm{m}\,\mathrm{s}^{-1}$ gives $\Delta x/c$ of 2.2 hours and, hence, a time
step of two hours will be chosen. If mRM is forced with hourly input, it
will aggregate the input over two consecutive time steps prior to the
routing. The calculated streamflow is then distributed to the previous
two time steps because these represent the mean flow over this
period. If aTS is forced with daily input, it will use internally 12
iterations of 2-hour time steps to route the water through the river
network. In this case, aTS will also return the average of the
calculated 12 streamflow values at each reach.

<a name="rMC"></a>

## Regionalized Muskingum-Cunge parametrization

The regionalized Muskingum-Cunge (rMC) parametrization implemented in
the mesoscale Hydrologic Model mHM calculates the Muskingum coefficients
$C_1$, $C_2$, and $C_3$ in the [routing recurrence](#mrm_recurrence) as a function of
high-resolution river network properties. The coefficients $C_1$, $C_2$,
and $C_3$ are parametrized as follows

$$
\begin{equation}
     C_1 = \nu_2; C_2 = \nu_1 - \nu_2; C_3 = 1 - \nu_1,
\end{equation}
$$

 where the parameters $\nu_1$ and $\nu_2$ are given as

$$
\begin{align}
     \nu_1 &= \frac{\Delta t}{\beta (1 - \epsilon) + \frac{\Delta t}{2}};\\
     \nu_2 &= \frac{\frac{\Delta t}{2} - \beta \epsilon}{\beta (1 - \epsilon) + \frac{\Delta t}{2}} \nonumber
\end{align}
$$

 following the nomenclature of appendix A2 in
 \cite SKA2010. This formulation is identical to
 the [coefficient expressions](#mrm_coefficients) given above, using
 $\beta = \Delta x / c$ in the equation above and substituting this
 equation into two equations above. The parameters
 $\beta$ and $\epsilon$ are then conceptualized as

$$
\begin{align}
     \beta &= \gamma_1 + \gamma_2 L + \gamma_3 S + \gamma_4 C ;\\
     \epsilon &= \gamma_5 \frac{S}{max(S)} , \nonumber
\end{align}
$$

 where $L$ is the length of the reach, $S$ is the slope of the reach,
 and $C$ is the fraction of impervious land cover within the floodplains
 (see table 4 in \cite KS2012). Overall, there are five global
 parameters $\gamma_1$ to $\gamma_5$ in the above equation that
 can be chosen by the user. The integration time step is fixed at one
 hour. To guarantee the numerical stability of the parameterization, the
 following upper and lower bounds are applied

$$
\begin{align}
   0 & < \epsilon \le 0.5,\\
   \frac{\Delta t}{2(1-\epsilon)} & < \beta \le \frac{\Delta t}{2 \epsilon},
\end{align}
$$

 where $\Delta t$ is set to one hour.
