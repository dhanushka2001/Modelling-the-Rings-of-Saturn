# Modelling the Rings of Saturn

This was my Final Year Project at the University of Warwick in 2022, which was given a [First-class](https://drive.google.com/file/d/1U0OKGjXay2MrkH1fxmdhtpEDRbPdHB9t/view?usp=sharing).

Degree: BSc Mathematics and Physics

[My thesis paper](https://drive.google.com/file/d/1HvYnVZAVaakcHQqrsdpO8QMVmfxp50X4/view?usp=sharing)

[PX319 Project List 2021-22](https://drive.google.com/file/d/1BXwNT60hRM-PJHbrFmTNUYFMOWn1Kz01/view?usp=sharing)

<!--
> more text

Hi $`m_s = 2`$

```math
SE = \frac{\sigma}{\sqrt{n}}
```
$$
m_{COM} = 2
$$

```math
\begin{aligned}
\dot{x} & = \sigma(y-x) \\
\dot{y} & = \rho x - y - xz \\
\dot{z} & = -\beta z + xy
\end{aligned}
```

Hi $`m_s = 2`$

<details><summary> Test </summary>

## Test
[<img width="2766" height="1364" alt="image" src="https://github.com/user-attachments/assets/55930fb5-a2d2-4aac-a539-142423533772" />](https://ciclops.org/view/2230/In-Saturns-Shadow---the-Pale-Blue-Dot.html)
_Figure 1: A picture of Saturn and its rings backlit by the Sun, taken by Cassini in 2006 (NASA/JPL-Caltech/Space
Science Institute)._
dsdsdsd

dsdsds Hi $`m_s = 2

```math
m_s = 2
```

```math
\begin{aligned}
\dot{x} & = \sigma(y-x) \\
\dot{y} & = \rho x - y - xz \\
\dot{z} & = -\beta z + xy
\end{aligned}
```

sasasas

sasasas

</details>
-->

## Authors
Dhanushka Jayagoda - <dhanushka2001@gmail.com>
### Affiliation
![Warwick University logo](images/Shield_of_the_University_of_Warwick-small.png)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;![Warwick University logo](images/WarwickLogo-small.png)   

<!--
## Section 1
01. text
   ![Saturn's ring gaps annotated](images/saturn4.png)
-->

## Final Report

<details><summary> Synopsis </summary>

## Synopsis

* The main objective of this project is to model the rings of Saturn by creating a
computer program, to see and explain how gaps form, as well as other structures
like spiral density waves, and wakes from embedded moons.

* The results will be
compared with observational data and theory to see how well they match.

* To model ring particles, the velocity Verlet method will be used to integrate Newton’s
equations of motion. Newton’s law of gravitation and Kepler’s third law will be
used.

* For satellite wakes, the ring particles will be modelled as an incompressible
fluid in a Keplerian shear flow, using 2D diffusion-advection equations. Gauss-Seidel iteration and the explicit Euler method will be used to solve the unsteady 2D
diffusion and advection equations.

* The results from the simulations show that gaps
form as predicted by the theory, in locations where ring particles are in resonance
with perturbing moons. From observational data, most gaps are further out than
what theory predicts, due to factors like the eccentricity of the moon’s orbit and
co-orbital moons which the simulations took into account. Results from the fluid
simulation show embedded moons do cause wakes.

</details>

<details><summary> Introduction </summary>

## Introduction
[<img width="2766" height="1364" alt="image" src="https://github.com/user-attachments/assets/55930fb5-a2d2-4aac-a539-142423533772" />](https://ciclops.org/view/2230/In-Saturns-Shadow---the-Pale-Blue-Dot.html)
_Figure 1: A picture of Saturn and its rings backlit by the Sun, taken by Cassini in 2006 (NASA/JPL-Caltech/Space
Science Institute)._


Saturn's rings are a majestic sight, with its radial symmetry, they exhibit many interesting features, from spiral density waves to the multitude of gaps. Some of these gaps are due to orbital resonances from Saturn's moons, one of the most notable being the Huygens gap, caused by the 2:1 inner Lindblad resonance (ILR) from Mimas.[^Paper_4] A ring particle in this location would have an orbital period half that of Mimas. Other gaps are caused by moons inside the rings, which produces a wake at the gap boundary. However, there are some gaps that are yet to be fully explained, some have proposed asteroid impacts,[^Paper_3] as well as interparticle collisions and excitation of density waves.[^Paper_2] This project aims to simulate these mechanisms, particularly, resonances, embedded moons, and spiral density waves, to see if they match with observation and theory.

[^Paper_4]: [K. Baillié et al., Formation of the Cassini Division – I. Shaping the rings by Mimas inward migration, Monthly Notices of the Royal Astronomical Society, 2933 (2019).](https://academic.oup.com/mnras/article/486/2/2933/5364573)
[^Paper_2]: [J. J. Lissauer, J. N. Cuzzi, Resonances in Saturn's Rings, Astronomical Journal **87**, 1052 (1982).](https://adsabs.harvard.edu/pdf/1982AJ.....87.1051L)
[^Paper_3]: [C. Kendrick, B. Kendrick, Computer Simulation of Saturn’s Ring Structure, 15 (2013).](https://drive.google.com/file/d/1jV_VsJOvCTXmZS-OiWrIcpOdBpUnkLyv/view?usp=sharing)

[<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/5ecf5417-5057-453a-b628-a17fd3dcbb67" />](https://ciclops.org/view/8489/The-Realm-of-Daphnis.html)
_Figure 2: A picture of Daphnis, a shepherding moon, first discovered by Cassini in 2005, taken by Cassini in 2017, leaving a wavy wake behind and in front of it (NASA/JPL-Caltech/Space Science Institute)._


Saturn's rings were first discovered in 1610 by Galileo, though he mistook them for two large satellites.[^Paper_1] Huygens was the first to propose it was a ring in 1655.[^Paper_1] In 1675, Cassini had discovered a gap in the ring,[^Paper_1] which we now call the "Cassini Division", he was also the first to propose that the rings were made out of smaller particles.[^Paper_1] In 1785, Laplace proved mathematically that even the narrow, rotating rings he envisioned would be unstable,[^Paper_1] however the idea that the rings were solid remained for almost another century. In 1857, Maxwell showed that the ring system had to be composed of many smaller particles, ending the 200-year-old belief that the rings were solid.[^Paper_1]

[^Paper_1]: [A. F. O'Donel, The Planet Saturn, 85, 92, 114, 119, 123, 187 (1962).](https://archive.org/details/planetsaturnhist0000alex/page/n9/mode/2up)

In recent years, NASA's Cassini probe has given us clear images of Saturn and its rings, which has led to the discovery of new moons, including shepherding moons which are embedded in the rings themselves, and leave a wake in their path. This has allowed many to study the gaps in great detail, and is what will be used to compare with the results.

Recent discoveries have found many exoplanets with ring systems much larger than Saturn's.[^Paper_5][^Paper_7] Understanding how Saturn's rings form can help to understand how much larger ring systems can form. Furthermore, the underlying physics behind ring structure formation, in particular, perturbation and resonance, is identical to that which forms spiral galaxies,[^Paper_6] as was first proposed by Goldreich and Tremaine,[^Paper_12] as well as accretion disk structure.[^Paper_8]

[^Paper_5]: [S. Rieder, M. A. Kenworthy, Constraints on the size and dynamics of the J1407b ring system, A&A **596**, A9 (2016).](https://www.aanda.org/articles/aa/full_html/2016/12/aa29567-16/aa29567-16.html)
[^Paper_6]: [T. C. Junqueira et al., A new model for gravitational potential perturbations in disks of spiral galaxies. An application to our Galaxy, A&A **550**, A91 (2013).](https://www.aanda.org/articles/aa/full_html/2013/02/aa19769-12/aa19769-12.html)
[^Paper_7]: [N. C. Santos et al., Detecting ring systems around exoplanets using high resolution spectroscopy: the case of 51Pegb, Earth and Planetary Astrophysics (2015).](https://www.aanda.org/articles/aa/full_html/2015/11/aa26673-15/aa26673-15.html)
[^Paper_8]: [O. Donmez, Numerical simulation of small perturbation on an accretion disk due to the collision of a star with the disk near the black hole, Astrophysics and Space Science **305**, 187 (2006).](https://arxiv.org/abs/astro-ph/0604249)
[^Paper_12]: [P. Goldreich, S. Tremaine, The excitation and evolution of density waves, Astrophysical Journal **222**, 850-858 (1978).](https://adsabs.harvard.edu/pdf/1978ApJ...222..850G)

In most cases, the perturbations from a moon do not accumulate on a particle, but get overpowered by the gravity of Saturn. However, if a particle and moon's orbital frequency are an integer ratio, the moon will periodically appear in the same location for the particle, which leads to gravitational perturbations accumulating, and causing the particle's orbit to become eccentric. This is one of the mechanisms that causes a visible gap to form in the ring.

Consider a 3-body system, with a primary body (Saturn), with mass $`m_p`$; an orbiting secondary body (moon), with mass $`m_s`$; and a ring particle. The particle will be orbiting the primary, while being perturbed by the secondary. The primary's frame is non-inertial, as it is orbiting around the system's centre-of-mass (COM). By definition,

<p align="center">
$\vec{r'_{COM}} = \frac{m_p\vec{r'_p}+m_s\vec{r'_s}}{m_p+m_s} = \vec{0},$ &nbsp; (1) <br>
$\implies \vec{r'_p} = -\frac{m_s}{m_p}\vec{r'_s},$ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (2) <br>
</p>

where $`\vec{r'_{COM}}`$, $`\vec{r'_p}`$ and $`\vec{r'_s}`$ are the positions of the COM, primary and secondary relative to the COM. $`m_p >> m_s \implies \vec{r'_p} \simeq \vec{0},`$ $`\therefore \vec{r'_p} \simeq \vec{r'_{COM}}`$. In reality, the mass of Saturn is $`\sim10^3`$ times heavier than its most massive moon, Titan, and $`\sim10^7`$ times heavier than the most massive moon used in this project, Mimas. For this reason, and for simplifying calculations, the primary will be approximated as an inertial frame.

The angle between $`\vec{r}`$ and $`\vec{r_s}`$ is $`\varphi=\theta - \theta_s`$, where $`\theta`$ and $`\theta_s`$ denote the particle and secondary's angle from the $`\vec{x}`$-axis. Switching to polar coordinates, $`\Phi_s`$ is a function of $`r`$ and $`\varphi`$, and crucially,  $`\Phi_s(r,\varphi)=\Phi_s(r,\textminus\varphi)`$. Because $`\Phi_s(r,\varphi)`$ is symmetric in $`\varphi`$, it can be Fourier expanded, with all the odd parts equal to zero, giving,[^Paper_9]

[^Paper_9]: [J. M. Hahn, The Dynamics of Planetary Systems and Astrophysical Disks, 112, 99 (2009).](https://gemelli.spacescience.org/~hahnjm/book/chap6.pdf)

<p align="center">
$\Phi_s(r,\varphi) = \frac{1}{2}\phi_0(r) + \sum_{m=1}^{\infty} \phi_m(r)\cos(m\varphi).$ &nbsp;&nbsp;&nbsp;&nbsp; (3) <br>
</p>

This Fourier expansion has infinite terms, each corresponding to a Lindblad resonance, which is a location where the cumulative perturbations from the secondary is large.[^Paper_9] For low $`m^{\textup{th}}`$-order terms, these locations are far apart, which means one term is dominating, and we can ignore the infinite other terms.[^Paper_9]

[<img width="4510" height="1794" alt="lindblad4" src="https://github.com/user-attachments/assets/b913fce7-97ec-44a5-8d7b-cba6731bd278" />](https://commons.wikimedia.org/wiki/File:Lindblad_resonance_sites.png)
_Figure 3: The radii of the secondary's $`m^{\textup{th}}`$ inner (ILR) and outer (OLR) Lindblad resonances relative to the secondary's orbit, located at the corotation circle (CC) (Dhanushka Jayagoda, 2022)._

Assuming that the secondary's orbit is circular, and the primary's gravitational potential is Keplerian, the locations of the Lindblad resonances are given by,[^Paper_9]

<p align="center">
$r = \left(1 - \frac{\epsilon}{m}\right)^\frac{2}{3}a_s = \left(\frac{m-\epsilon}{m}\right)^\frac{2}{3}a_s,$ &nbsp;&nbsp;&nbsp;&nbsp; (4) <br>
</p>

where $`\epsilon=\pm1`$, and $`a_s`$ is the secondary's orbital radius. Fig. 3 shows the locations of the Lindblad resonances, $`r=[(m-1)/m]^{2/3}a_s`$ are the ILRs, and $`r=[(m+1)/m]^{2/3}a_s`$ are the OLRs. Because of the assumptions, the locations will slightly differ if the secondary's orbit is eccentric, or if the primary body is oblate.[^Paper_9]

From Kepler's third law we have,

<p align="center">
$\left(\frac{T}{T_s}\right)^2 = \left(\frac{R}{R_s}\right)^3,$ &nbsp;&nbsp;&nbsp;&nbsp; (5) <br>
</p>

where $`T`$, $`T_s`$, $`R`$ and $`R_s`$ denote the particle and secondary's orbital periods and radii. Rearranging and solving for $R$ gives,

<p align="center">
$R = \left(\frac{T}{T_s}\right)^\frac{2}{3}R_s.$ &nbsp;&nbsp;&nbsp;&nbsp; (6) <br>
</p>

Letting $`k=T/T_s`$ denote the ratio of the orbital periods, we can see that Eq. (4) and \(6) are similar. From Eq. (4) and (6), we have that the ratio of the orbital periods for a $`m^{\textup{th}}`$-order Lindblad resonance is $`k=(m-\epsilon)/m`$. For a $`m=2`$ ILR, $`k=\frac{1}{2}`$, meaning the resonance location occurs where a particle has an orbital period half that of the secondary, as expected. Since the majority of ring particles are located on the inside of Saturn's moons' orbits, only ILRs will be considered in this project.
 
Saturn also features two moons, Janus and Epimetheus, that are in a co-orbital configuration, which is the only configuration of its kind with two moons, in the Solar System.[^Paper_15] Janus and Epimetheus have orbital radii that are closer than their own diameters, yet they do not collide, but rather periodically exchange momentum and switch orbital radii every 4 years.[^Paper_15] Janus, the more massive of the two, is said to create the most visible spiral density waves from its excited 2:1 Lindblad resonance site, due to its oscillatory orbit.[^Paper_16]

[^Paper_15]: [M. El Moutamid et al., How Janus’ orbital swap affects the edge of Saturn’s A ring?, Icarus **279**, 125, (2016).](https://arxiv.org/abs/1510.00434)
[^Paper_16]: [J. E. Colwell et al., Saturn from Cassini-Huygens, **13**, 392, (2009).](https://link.springer.com/book/10.1007/978-1-4020-9217-6)

</details>

<details><summary>Computational details</summary>

## Computational details

### A. Numerical methods: Euler and Verlet

The two methods considered for integrating Newton's equations of motion in this project are, Euler's method and the velocity Verlet method.

Euler's method is as follows,

1. Calculate $`\vec{r}(t+\Delta t) = \vec{r}(t) + \vec{v}(t)\Delta t + \frac12\vec{a}(t)\Delta t^2,`$
2. Work out $`\vec{a}(t+\Delta t)`$ from the gravitational interactions using $`\vec{r}(t+\Delta t),`$
3. Calculate $`\vec{v}(t+\Delta t) = \vec{v}(t) + \vec{a}(t)\Delta t.`$

The Velocity Verlet algorithm is as follows,

1. Calculate $`\vec{r}(t+\Delta t) = \vec{r}(t) + \vec{v}(t)\Delta t + \frac12\vec{a}(t)\Delta t^2,`$
2. Work out $`\vec{a}(t+\Delta t)`$ from the gravitational interactions using $`\vec{r}(t+\Delta t),`$
3. Calculate $`\vec{v}(t+\Delta t) = \vec{v}(t) + \frac12\left(\vec{a}(t) + \vec{a}(t+\Delta t)\right)\Delta t.`$

The only difference between the methods is the final step, where the Verlet method combines the old and new acceleration vectors, instead of just the old acceleration vector. This change drops the magnitude of $`\%error`$ by a factor of $`10^3`$ for a timestep $`\Delta t=1000~s`$, and changes the order of error from $`\mathcal{O}(\Delta t)`$ to $`\mathcal{O}(\Delta t^2)`$, as can be seen from Fig. 4.

<img width="1920" height="977" alt="saturn_error2 1" src="https://github.com/user-attachments/assets/924e9695-b781-424a-b619-26d44bb91ed1" />
<img width="1920" height="977" alt="saturn_error6 1" src="https://github.com/user-attachments/assets/7d38e504-7791-4af6-b75f-cd0eb055d87f" />

_Figure 4: Above is two graphs plotting the $`\%error`$ in position and velocity of Mimas against timestep ($`s`$) for both Euler and Verlet methods. Below is a graph of $`\%error`$ in position of Mimas against timestep using Verlet. Mimas is modelled to have a circular orbit. (Dhanushka Jayagoda, 2022)_
    
Because of the improvement in stability, velocity Verlet will be the method used in this project. For both methods, Newton's law of gravitation is used to calculate the gravitational interactions from Saturn and its moons, and Kepler's third law is used, to calculate orbital periods and initial velocities, given an initial orbital radius, $`R`$.

From Kepler's third law, we are given that,[^Paper_14]

[^Paper_14]: [R. R. Bate, D. D. Mueller, J. E. White, Fundamentals of astrodynamics, 9, 193, 221 (1971).](https://cmp.felk.cvut.cz/~kukelova/pajdla/Bate,%20Mueller,%20and%20White%20-%20Fundamentals%20of%20Astrodynamics.pdf)

<p align="center">
$T^2 = \left(\frac{4\pi^2}{Gm_p}\right)R^3,$ &nbsp;&nbsp;&nbsp;&nbsp; (7) <br><br>
$\therefore &nbsp;&nbsp;&nbsp;&nbsp; T = 2\pi\sqrt{\frac{R^3}{Gm_p}}.$ &nbsp;&nbsp;&nbsp;&nbsp; (8) <br>
</p>

Substituting Eq. (8) into the formula for tangential velocity gives,

<p align="center">
$v = \frac{2\pi R}{T} = \sqrt{\frac{Gm_p}{R}},$ &nbsp;&nbsp;&nbsp;&nbsp; (9) <br>
</p>

where $G$ is the gravitational constant. Acceleration is calculated using Newton's law of gravitation, for a ring particle $`i`$,[^Paper_14]

<p align="center">
$\ddot{\vec{r_i}} = G \cdot \sum_{j\neq i} \frac{m_j}{r_{ij}^3} \vec{r_{ij}},$ &nbsp;&nbsp;&nbsp;&nbsp; (10) <br>
</p>

With $`n`$ ring particles, the total interactions amounts to $`n(n-1)/2`$, which is of order $`\mathcal{O}(n^2)`$. In the interest of speed, the particles in the Verlet method simulation will be modelled as non-interacting, only interacting with Saturn and its moons, which is of order $`\mathcal{O}(n)`$ interactions. In reality, ring particles will experience self-gravity with neighbouring particles, however Keplerian shear mostly prevents this from becoming dominant.[^Paper_13]

[^Paper_13]: [C. Murray, Saturn's dynamical rings, Physics Today **60**, 8, 74 (2007).](https://physicstoday.aip.org/quick-study/saturns-dynamical-rings)

### B. Special perturbation methods: Cowell and Encke

<img width="3000" height="3000" alt="enckemethod" src="https://github.com/user-attachments/assets/c3eeb295-1e68-439e-a3ee-d63d43a5582c" />

_Figure 5: A diagram showing Encke's method, with $`n\_subcycles=3`$. Perturbations are exaggerated. (Dhanushka Jayagoda, 2022)_

There are two methods that can be used when calculating perturbations, Cowell's method and Encke's method. Cowell's method adds the perturbing interactions together and integrates in time continuously, whereas Encke's method involves subcycles, where perturbations are accumulated, and the unperturbed, osculating orbit is used as a reference, for which it rectifies discretely at every cycle. Fig. 5 shows how all the perturbations from the last rectification are accumulated, and added after the end of the current cycle. Perturbations of velocity will also be accumulated, as well as of position, due to the fact that $`\vec{r}(t+\Delta t)`$ and $`\vec{v}(t+\Delta t)`$ both need $`\vec{a}(t)`$. For the Verlet method, perturbations are calculated as follows,

<p align="center">
$d\vec{r} = \frac12 \Delta t^2 \sum_{n=0}^{N-1} \vec{a_s}(t+n\Delta t) = \sum_{n=1}^{N} d\vec{r_n},$ &nbsp;&nbsp;&nbsp;&nbsp; (11) <br><br>
$d\vec{v} = \frac12 \Delta t \sum_{n=0}^{N-1} \big(\vec{a_s}(t+n\Delta t) + \vec{a_s}(t+n\Delta t + \Delta t)\big) = \sum_{n=1}^{N} d\vec{v_n},$ &nbsp;&nbsp;&nbsp;&nbsp; (12) <br>
</p>

where $`N=n\_subcycles`$ is the number of subcycles, and $`t`$ is the simulation time of the current cycle. These perturbations are for the next cycle, at $`T=t+N\Delta t`$.
    
Cowell's method is better for when you want to find the path of a particle at a particular time, since the program updates continuously, however Encke's method allows you to go much further in time with the sacrifice that the data points are no longer continuous but occur at discrete jumps in time. Encke's method is also faster, however it is less accurate when a particle experiences large perturbations, such that $`dr/r`$ is no longer small.[^Paper_14]
    
The total simulation time,

<p align="center">
$t_s = \Delta t \cdot n\_subcycles \cdot n\_cycles,$ &nbsp;&nbsp;&nbsp;&nbsp; (13) <br>
</p>

where $`n\_cycles`$ is the number of cycles. Note that Cowell's method is a case of Encke's method with $`n\_subcycles=1`$. One of the known problems with Cowell and Encke's methods is the exponential error that accumulates over time. The only solutions are to reduce the timestep, or avoid particles getting too close to the moons.[^Paper_15]

### C. The elliptic orbit

For a moon with a circular orbit, radius $`r`$, its positions can simply be found through the parametric equation, $`(x,y)=\big(r\cos(2\pi t/T),r\sin(2\pi t/T)\big)`$, where $`t`$ and $`T`$ are the simulation time and orbital period. Modelling an elliptic orbit is not as simple. Parametric equations for an elliptic orbit can however, be approximated, by iteratively solving Kepler's equation.
    Kepler's equation is given by $`M = E - e\sin(E)`$, where $`M`$ is known as the "mean anomaly", $`e`$ is the eccentricity, and $`E`$ is the "eccentric anomaly".[^Paper_14] Calculating $`M`$ given $`E`$ is trivial, however the inverse problem is not, since the equation is transcendental, which means an approximate solution will need to be found instead.[^Paper_14] The steps for calculating an approximate position are as follows,[^Paper_14]

1. Calculate the mean motion, $`n = 2\pi/T`$, $`T=`$ orbital period,
2. Calculate the mean anomaly, $`M = nt`$, $`t=`$ time since perihelion (point on elliptic orbit closest to primary),
3. Compute the eccentric anomaly, $`E`$, by using Newton's method on Kepler's equation. Let,

   <p align="center">
   $f(E) = E - e\sin(E) - M(t) = 0,$ &nbsp;&nbsp;&nbsp;&nbsp; (14) <br>
   $E_{n+1} = E_{n} - \frac{f(E_n)}{f'(E_n)} = E_n - \frac{M_n-M(t)}{1-e \cos(E_n)},$ &nbsp;&nbsp;&nbsp;&nbsp; (15) <br>
   </p>

    where $`M_n=E_n-e\sin(E_n)`$.[^Paper_14] Eq. (15) can be iterated as many times as needed, until $`M_n - M(t)`$ becomes small enough.[^Paper_14] An initial value, $`E_0=\pi`$ is recommended for stability, however for $`e \approx 0`$, $`E_0=M(t)`$ is ample.[^Paper_14]

4. Calculate the approximate new position $`(x,y)=\left(a(\cos(E)-e),b\sin(E)\right)`$, where $`a`$ is the semi-major axis, and $`b=a\sqrt{1-e^2}`$ is the semi-minor axis.
   
Almost all of Saturn's moons can be modelled much faster this way, rather than using a numerical method. However, the one exception to this is the Janus-Epimetheus co-orbital configuration, which needs to be simulated with the Verlet method, to accurately model the horseshoe orbit.[^Paper_15]

### D. 2D diffusion-advection

In order to model the wake trails left by embedded moons, interparticle collisions can no longer be ignored. A small region near the moon will be considered, and area discretised so that a $`4~km`$-by-$`4~km`$ region is represented by one grid square, with a unique density, $`\rho_i`$, and velocity, $`\vec{v_i}`$. Gauss-Seidel iteration will be used to solve the Forward-Time Centered-Space (FTCS) approximation of the unsteady 2D diffusion equation.[^Paper_17] The gravitational interaction from the satellite, modelled in the centre of the grid, will evolve the velocity field and the explicit Euler method will be used to solve the advection equations. For the boundary conditions, density will be allowed to flow out of the system, as source terms will be included in the sides to account for the Keplerian shear flow. The frame of the moon is non-inertial, as it is orbiting and rotating, which keeps Saturn directly above the system. Given the region is small, the fictitious force from the frame orbiting will cancel with the centripetal force from the particles, and any horizontal component remaining is negligible compared to the interaction with the moon. Likewise for the fictitious forces arising from the frame rotating.

[^Paper_17]: [J. D. Hoffman, Numerical Methods for Engineers and Scientists, 546, (1992).](http://freeit.free.fr/Finite%20Element/Hoffman,_Numerical_Methods_for_Engineers&Scientists,2001.pdf)

</details>

<details><summary>Results</summary>

## Results

### A. Resonance strengths

<img width="5917" height="4530" alt="saturns_moons_strengths4" src="https://github.com/user-attachments/assets/7e05531a-08f4-44be-a28d-9a6adb086566" />

_Figure 6: Magnitude of the accumulated accelerations from each moon at varying distances from Saturn. (Dhanushka Jayagoda, 2022)_

<img width="410" alt="resonance_strengths_dt=1000_k=10000" src="https://github.com/user-attachments/assets/f16b2893-41ac-4cdc-b366-8eb766370958" />
<img width="410" alt="resonance_strengths_dt=1000_k=100000" src="https://github.com/user-attachments/assets/8194fc75-fa97-4ea2-a1ad-6dbb4b39aa3b" />

_Figure 7: Magnitude of the accumulated accelerations, with all the moons' perturbations combined, at varying distances from Saturn. Total simulation time: $`\sim0.32 ~yrs`$ (left) and $`\sim3.17 ~yrs`$ (right). (Dhanushka Jayagoda, 2022)_

In order to measure the strengths of the Lindblad resonances from each other moons, and get a rough idea of which moons matter the most, $`10,000`$ particles were placed between $`74,500~km`$ and $`140,220~km`$ from Saturn, covering Saturn's A, B, C, and F rings, along with 14 of Saturn's most prominent moons, each interacting with all $`10,000`$ particles. The particles and moons were simulated with fixed circular orbits, and the perturbations from each of the moons are separately accumulated. Since the particles were fixed in a circular orbit, no numerical method was needed to calculate the new positions. The simulation ran using Cowell's method, with $`10,000`$ cycles, each with a timestep $`\Delta t=1,000~s`$, amounting to $`\sim0.32~yrs`$ simulation time. Fig. 6 shows the results from this; by far the largest spike is from the 2:1 Mimas ILR, followed by the Janus-Epimetheus ILRs. Some of the embedded moons like Prometheus, Pan and Atlas show large perturbations, however these can't be considered as single resonance sites, as the particles close to these moons feel a contribution of many higher $`m^{th}`$-order ILRs. It is interesting to note that Titan, Saturn's heaviest moon, which is $`\sim10^4`$ times heavier than Mimas, does not even appear on the graph, this is due to its lowest order ILR, the 2:1 resonance, being located at $`r=(1/2)^{2/3}a_s \approx 769,730~km`$, far clear of the rings. From Fig. 7 it can be seen that over time, regions not close to any Lindblad resonances dampen in perturbations, whereas sites like the 2:1 Janus-Epimetheus ILR, spike in accumulated acceleration after more time has passed.

### B. Huygens gap

<img width="410" alt="saturnringdensity6graph2" src="https://github.com/user-attachments/assets/c465363a-a891-444b-ab42-e94fc4c7e84b" />
<img width="410" alt="saturnringdensity6graph2_ellipse" src="https://github.com/user-attachments/assets/f1efdab3-6ee2-4e50-9537-c132116e3423" />

_Figure 8: Ring particle density over time. $`10,000`$ particles initially between $`115,500~km`$ and $`118,000~km`$, with just Mimas and Saturn. Mimas' orbit modelled as a circle (left) and an ellipse (right). Total simulation time: $`\sim25.8~yrs`$. (Dhanushka Jayagoda, 2022)_

From the investigation into resonance strengths, it is clear that the 2:1 Mimas ILR is the most dominant resonance affecting the rings. Mimas' semi-major axis, $`a_s=185,539~km`$, therefore its 2:1 ILR location should be at $`r\approx116,882~km`$. Placing $`10,000`$ particles between $`115,500~km`$ and $`118,000~km`$, with just Mimas and Saturn influencing them, and using the Verlet and Encke methods to orbit the particles, produces Fig. 8. The top plot was with Mimas' modelled with a circular orbit, and the bottom plot was Mimas' modelled with an elliptic orbit, accurate to its true eccentricity, $`e=0.0196`$. Fig. 8 shows how the gaps get slightly pushed out with an eccentric Mimas orbit, Fig. 10 shows how the radial velocity spikes further out with an eccentric orbit. In the circular orbit plot, there appears to be a ringlet inside the gap, however with the elliptic orbit, there is a clean gap as predicted by the theory. In reality however, the Huygens gap is located at a radial distance, $`r=117,580~km`$, as seen in Fig. 9. This anomaly could be explained by the oblateness of Saturn, which is said to shift Lindblad resonances outward by,[^Paper_9]

<p align="center">
$\Delta r_r \simeq \frac{J_2}{2}\left(\frac{m+\epsilon}{m-\epsilon}\right)\frac{R_p^2}{r_r},$ &nbsp;&nbsp;&nbsp;&nbsp; (16) <br>
</p>

where $`J_2`$ is the primary's zonal harmonic and $`R_p`$ is the equatorial radius of the oblate planet. Both the circular and elliptic plots appear to show density waves initially, at around $`116,000~km`$, and start to drift away from the resonance site as predicted by theory.[^Paper_9]

[<img width="410" alt="image" src="https://github.com/user-attachments/assets/44ce3c68-68ff-40de-9e01-9bce62cadac0" />](https://ciclops.org/view/2127/The-Huygens-Gap.html)
[<img width="410" alt="image" src="https://github.com/user-attachments/assets/6277d0cf-1afc-4e38-b300-f2981c77d512" />](https://ciclops.org/view/3858/Expanse-of-Ice.html)
_Figure 9: Image of the Huygens gap taken by Cassini in 2006 (left), annotated image of the Cassini division with the Huygens gap, taken in 2007 by Cassini (right) (NASA/JPL-Caltech/Space Science Institute)._

<img width="1919" height="910" alt="radial_v" src="https://github.com/user-attachments/assets/aa0fcc80-d485-4077-b9a9-68fd38f05377" />

_Figure 10: The two diagrams above show the difference in radial velocity for eccentricity $`e=0`$ (left), and $`e=0.0196`$ (right). Simulation time: $`\sim25.8~yrs`$. (Dhanushka Jayagoda, 2022)_

### C. Co-orbital moons

[<img width="3849" alt="radial_v" src="https://github.com/user-attachments/assets/9643c067-8f26-4752-ab0f-ec3e10f16908" />](https://www.youtube.com/watch?v=P12ifsFwC6k)
_Figure 11: Locations of Janus and Epimetheus in a top-down view. The diagram above shows the horseshoe orbit when viewed in a rotating frame, and also shows them at their closest meeting, $`\sim15,000~km`$. (Click to watch animation). (Dhanushka Jayagoda, 2022)_

Janus and Epimetheus are two moons that orbit in roughly the same radius, swapping orbits and exchanging momentum periodically every 4 years. In order to model this, the Verlet method had to be used, as well as the Cowell method to get every position continuously rather than after discrete intervals in time. Initialising Janus to $`r=151,460~km`$, and Epimetheus to $`r=151,410~km`$ produced a 4-year periodicity, matching observational data.[^Paper_15]

<img width="410" alt="januslinegraph5" src="https://github.com/user-attachments/assets/231eb220-1e20-4ab9-a651-6f49bcc22a8c" />
<img width="410" alt="saturnringdensitygraph_janus_epimetheus2" src="https://github.com/user-attachments/assets/e52409f9-940f-4bc7-8ebe-970551460dd8" />

_Figure 12: Line graphs of radial positions (left) and ring particle density (right) for 10,000 particles close to the 2:1 ILR of Janus. Total simulation time: $\sim25.8$yrs. (Dhanushka Jayagoda, 2022)_

[<img width="1024" alt="image" src="https://github.com/user-attachments/assets/2f399781-edbb-47a1-9cb5-6ee184b0903a" />](https://ciclops.org/view/8612/Staggering-Structure.html)
_Figure 13: A spiral density wave at the 2:1 Janus ILR (NASA/JPL-Caltech/Space Science Institute)._

With the horseshoe orbit modelled, it is possible to influence the $`10,000`$ particles realistically. The top graph of Fig. 12 shows how particles start to clump together initially, and there's a clear drop in density at the radial distance $`r=95,900~km`$. The expected 2:1 ILR location for Janus from the theory is $`\sim95,414~km`$, this outward movement of the resonance location is indicative of the unique orbit pattern of Janus and Epimetheus. The bottom graph shows clearly how the 2:1 line is below the actual resonance location. The top graph of Fig. 12 appears to exhibit density waves as shown in Fig. 13, which seem to decay as time progresses.

### D. Shepherding moon

<img width="1801" height="613" alt="foo000" src="https://github.com/user-attachments/assets/176e6afd-a273-4dab-9c54-24e440079e50" />
<img width="1801" height="613" alt="foo020" src="https://github.com/user-attachments/assets/97023021-a57a-411a-a352-67fc566eca05" />
<img width="1801" height="613" alt="foo040" src="https://github.com/user-attachments/assets/e06e10c5-d34a-4af1-b76f-9141d88c3df6" />
<img width="1801" height="613" alt="foo060" src="https://github.com/user-attachments/assets/06e021c1-666c-435e-afa0-d98f1831d665" />
<img width="1801" height="613" alt="foo080" src="https://github.com/user-attachments/assets/9eeba881-eaaa-4811-b514-e2a99214062d" />
<img width="1801" height="613" alt="foo100" src="https://github.com/user-attachments/assets/c3d7ae90-38ca-4305-9484-88f9b2b1d606" />
<img width="1801" height="613" alt="foo200" src="https://github.com/user-attachments/assets/6fbeda26-6e8e-4031-8527-873cdafe3592" />
<img width="1801" height="613" alt="foo300" src="https://github.com/user-attachments/assets/0b8cb1f2-5423-4363-90d5-603474b11d9b" />

_Figure 14: 2D diffusion-advection simulation for Daphnis over time. (Dhanushka Jayagoda, 2022)_

Daphnis is the shepherding moon investigated in this project. From Fig. 14 wave-like wakes are visibly created by the embedded moon. There is slight deviation to what appears in real photos like Fig. 2 and movies taken by Cassini, as the wake seems to appear stationary, however the wake produced in the simulation seems to be periodically generated by a vortex of ring particle density hitting the gap boundaries.

</details>

<details><summary>Conclusions</summary>

## Conclusions

In conclusion, the work done in this project has shown how Lindblad resonances and perturbation can create visible gaps in the rings, and how elliptic orbits, and co-orbital moon configurations can lead to gaps moving further out, closer to the perturbing moon. This project has also shown that there are still some discrepancies with observation, meaning there is more at work causing the gaps to move further from the theory predicted locations. To take this project further, the oblateness of Saturn would have to be investigated. Precession of moon orbits could also be considered, as well as drag caused from Saturn's atmosphere, which would create the faint D-ring which the current simulation doesn't show. The mechanisms behind moonlet formation would be something to look into, self-gravity between neighbouring ring particles can lead to spiral density waves as well as bending waves in the vertical direction, which was not touched upon in this project.

</details>

## License
GNU General Public License v3.0

## Citations
