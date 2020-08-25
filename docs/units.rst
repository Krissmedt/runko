Runko Unit System and Particle Integration
##########################################

Introduction
============
For an arbitrary system of units, the Maxwell equations can be written,

.. math::
\nabla\cdot\mathbf{E}=4\pi k_{1}\rho,


.. math::
\frac{\partial\mathbf{E}}{\partial t}=\frac{c{{}^2}}{k_{2}}\nabla\times\mathbf{B}-4\pi k_{1}\mathbf{J},


.. math::
\nabla\cdot\mathbf{B}=0,


.. math::
\frac{\partial\mathbf{B}}{\partial t}=-k_{2}\nabla\times\mathbf{E},


where :math:`k_{1}` and :math:`k_{2}` are constants corresponding to a given
unit system (cite jackson). For CGS (Gaussian) units, the constants
default to :math:`k_{1}=1` and :math:`k_{2}=c`

.. math::
\nabla\cdot\mathbf{E}=4\pi\rho,


.. math::
\frac{\partial\mathbf{E}}{\partial t}=c\nabla\times\mathbf{B}-4\pi\mathbf{J},


.. math::
\nabla\cdot\mathbf{B}=0,


.. math::
\frac{\partial\mathbf{B}}{\partial t}=-c\nabla\times\mathbf{E},


with the corresponding Lorentz acceleration

.. math::
\frac{d\mathbf{u}}{dt}=\frac{q}{m}\mathbf{E}+\frac{\mathbf{v}}{c}\times\mathbf{B}.


The unit system in runko follows the convention used in the PIC code
TRISTAN-MP (Buneman 1993). 
The variables are normalised w.r.t. appropriate
fiducial values. The most noteworthy scalings are those for time and
distance, which are normalised w.r.t. the grid spacing :math:`\mathbf{x}=\hat{\mathbf{x}}\Delta x`
and time-step :math:`t=\hat{t}\Delta t`. As velocity is the derivative
of time with distance, it follows that the system velocity is expressed
in terms of the spatial and temporal step sizes. The coordinate/proper
velocity is given by

.. math::
\mathbf{v}=\frac{d\mathbf{x}}{dt}=\frac{d\hat{\mathbf{x}}}{d\hat{t}}\frac{\Delta x}{\Delta t}=\hat{\mathbf{v}}\frac{\Delta x}{\Delta t},


with the corresponding four-velocity

.. math::
\mathbf{u}=\hat{\mathbf{u}}\frac{\Delta x}{\Delta t}=\hat{\mathbf{v}}\gamma(\mathbf{v})\frac{\Delta x}{\Delta t}.


To maintain stability and ensure the fidelity of the simulation, the
fraction :math:`\Delta x/\Delta t` is set to equal :math:`c` normalised as before
so that :math:`\Delta x/\Delta t=c/\hat{c}`. In this manner, effective
light speed can be set in the simulation via :math:`\hat{c}` and the CFL
condition strictly enforced at all times.
In all the formulas solved by the computer the :math:`\hat{c}` ends up playing the role of a time-step because it effectively reduces maximum signal velocity (i.e., speed of light) so that :math:`\Delta t = c \Delta x \hat{c}`.
Setting grid size and light speed to unity, :math:`\Delta x = 1` and :math:`c = 1`, the definition reduces to :math:`\hat{c} = \Delta t`.

The fields are scaled via the fiducial value :math:`B_{0}`, so :math:`\mathbf{E}=\hat{\mathbf{E}}B_{0}`
and :math:`\mathbf{B}=\hat{\mathbf{B}}B_{0}`. Charge, mass and current
are similarly scaled, giving :math:`q=\hat{q}q_{0}`, :math:`m=\hat{m}m_{0}`
and :math:`\mathbf{J}=\hat{\mathbf{J}}J_{0}`where all numerical variables
are denoted with a hat.

Applying these conventions to the evolution equations, which are needed
to advance the PIC system, the runko update equations take the following
form in code units

.. math::
\Delta[\hat{\mathbf{E}}]_{t}=c\frac{\Delta t}{\Delta x}\Delta[\hat{\mathbf{B}}]_{x}-\frac{4\pi J_{0}}{B_{0}}\Delta t\hat{\mathbf{J}}=\hat{c}\Delta[\hat{\mathbf{B}}]_{x}-\hat{\mathbf{J}},


.. math::
\Delta[\hat{\mathbf{B}}]_{t}=-c\frac{\Delta t}{\Delta x}\Delta[\hat{\mathbf{E}}]_{x}=-\hat{c}\Delta[\hat{\mathbf{E}}]_{x},


.. math::
\Delta[\hat{\mathbf{u}}]_{t}=\frac{\hat{q}q_{0}\hat{c}}{\hat{m}m_{0}c}B_{0}\Delta t\left(\hat{\mathbf{E}}+\frac{\hat{\mathbf{v}}c}{c\hat{c}}\times\hat{\mathbf{B}}\right)=\frac{\hat{q}}{\hat{m}}\left(\hat{\mathbf{E}}+\frac{\hat{\mathbf{v}}}{\hat{c}}\times\hat{\mathbf{B}}\right),


Note here for the velocity update that :math:`\Delta[\hat{\mathbf{u}}]_{t}`
is a numerical expression for acceleration and is thus scaled using
`\Delta x/\Delta t^{2}`, which is changed to :math:`\Delta tc/\hat{c}`
and multiplied to the right hand side. The update equation for position
is defined by the integral of proper velocity :math:`\hat{\mathbf{v}}`, which takes the following form once scalings are introduced

.. math::
\Delta[\hat{\mathbf{x}}]_{t}=\frac{\Delta t}{\Delta x}\frac{\hat{\mathbf{u}}}{\gamma(\hat{\mathbf{u}})}\frac{\Delta x}{\Delta t}=\frac{\hat{\mathbf{u}}}{\gamma(\hat{\mathbf{u}})}.


In all equations above, :math:`\Delta[\cdot]_{t}` and :math:`\Delta[\cdot]_{x}` indicate
an appropriate finite difference operator of form :math:`x_{n+1}-x_{n}` or :math:`t_{n+1/2}-t_{n-1/2}`. Note the missing step sizing due to the normalisation of the system.
To simplify the equations, the following ratios between the different fiducial values were defined:

.. math::
B_{0}=\frac{m_{0}c}{q_{0}\hat{c}\Delta t}=\frac{m_{0}c}{q_{0}\hat{c}}\left(\frac{c}{\hat{c}\Delta x}\right),


.. math::
\frac{J_{0}}{B_{0}}=\frac{1}{4\pi\Delta t},

which effectively make eqs. (\ref{eq:ampere_updateE}) - (\ref{eq:lorentz_update}) unitless. The current is then scaled as

.. math::
\mathbf{J}=J_{0}\hat{\mathbf{J}}=\hat{q}\frac{q_{0}}{\Delta x^{3}}\hat{\mathbf{v}}\frac{\Delta x}{\Delta t}=\hat{q}\hat{\mathbf{v}}\frac{q_{0}}{\Delta x^{2}\Delta t},

which implies :math:`J_{0}=q_{0}/\Delta x^{2}\Delta t`. The formula for :math:`J` can be thus be understood as a charge flux through a surface area of size :math:`\Delta x^2`. Substituting in the values for :math:`J_{0}=q_{0}/\Delta x^{2}\Delta t` yields

.. math::
\frac{q_{0}}{\Delta x^{2}}=\frac{m_{0}c}{4\pi q_{0}\hat{c}\Delta t},


which can be rearranged for :math:`\Delta x`to form

.. math::
\Delta x=4\pi\frac{q_{0}^{2}}{m_{0}}\left(\frac{\hat{c}}{c}\right)^{2},


using the substitution :math:`\Delta x/\Delta t=c/\hat{c}.`

The particles simulated in PIC generally do not correspond to physical particles. Instead they represent groups/clusters of physical particles traveling together with same velocity and direction. These "macroparticles" thus allows the model to simulate a plasma of a given mass, charge and size without 1:1 resolution of the physical particles.
For a macroparticle containing :math:`N` electrons/positrons, the reference charge will equal :math:`q_{0}=Nq_{e}` and the expression becomes

.. math::
\Delta x=4\pi N\frac{q_{e}^{2}}{m_{e}}\left(\frac{\hat{c}}{c}\right)^{2}=4\pi N\hat{c}^{2}r_{e},


where :math:`m_{e}`, :math:`q_{e}=e` are the electron rest mass and elementary
charge respectively, resulting in the appearance of the classical
electron radius :math:`r_{e}=e^{2}/m_{e}c^{2}\approx2.82\cdot10^{-13}`.
The field scaling in Gaussian units then becomes

.. math::
B_{0}=\frac{m_{e}c}{e\hat{c}}\left(\frac{c}{\hat{c}\Delta x}\right)=\frac{e}{r_{e}\hat{c}^{2}\Delta x}.


Converting from code to Gaussian units now just requires selecting a grid length scale :math:`\Delta x`

.. math::
\mathbf{B}=\hat{\mathbf{B}}B_{0}=1.705\cdot10^{3}\frac{\hat{\mathbf{B}}}{\hat{c}^{2}}\,\frac{1\,\mathrm{cm}}{\Delta x}\mathrm{G},


.. math::
\mathbf{E}=\hat{\mathbf{E}}B_{0}=1.705\cdot10^{3}\frac{\hat{\mathbf{E}}}{\hat{c}^{2}}\,\frac{1\,\mathrm{cm}}{\Delta x}\mathrm{statvolt\,cm^{-1}},


.. math::
\mathbf{J}=\frac{B_{0}}{4\pi\Delta t}\hat{\mathbf{J}}=\frac{ec}{4\pi r_{e}\hat{c}{{}^3}\Delta x{{}^2}}\hat{\mathbf{J}}\approx4.056\cdot10^{12}\frac{\hat{\mathbf{J}}}{\hat{c}{{}^3}}\,\left(\frac{1\,\mathrm{cm}}{\Delta x}\right)^{2}\,\mathrm{statcoul\,s^{-1}},


.. math::
q=\hat{q}\frac{e\Delta x}{4\pi\hat{c}^{2}r_{e}}\approx1.356\cdot10^{2}\frac{\hat{q}}{\hat{c}^{2}}\,\frac{\Delta x}{1\,\mathrm{cm}}\,\mathrm{statcoul},

Using the above fundamental scalings, expressions for the conversion
of the derived plasma quantities from code units can also be determined.
The relativistic plasma frequency is given by

.. math::
\omega_{p}^{2}=\omega_{p,-}^{2}+\omega_{p,+}^{2}=\frac{4\pi q^{2}n}{<\gamma>m_{e}}\left(1+\frac{m_{-}}{m_{+}}\right),


.. math::
\hat{\omega}_{p}^{2}=\omega_{p}^{2}\Delta t^{2}=\frac{\hat{q}^{2}N_{ppc}}{\langle \gamma \langle \hat{m}}\left(1+\frac{m_{-}}{m_{+}}\right)=\frac{>|\hat{q}|N_{ppc}}{<\gamma>}\left(1+\frac{m_{-}}{m_{+}}\right),

for a 2-species plasma with :math:`+` and :math:`-` particles and masses :math:`m_-` and :math:`m_+`. In the case of electron-positron pair plasma :math:`m_- = m_+`, but for ions a mass scale :math:`m_+/m_-` is requiredthat in a real plasma is :math:`m_i/m_e`. The parameter :math:`N_{ppc}`denotes the number of particles per cell per species.

The numerical plasma frequency can then be set by selecting a skin depth resolution 

.. math::
R_p = \frac{c/\omega_p}{\Delta x},

which expresses the length of a plasma oscillation at velocity :math:´c´ non-dimensionalised w.r.t. the grid length.

.. math::
\hat{\omega}_{p,0}=\frac{\hat{c}\Delta x}{R_{p}\Delta x}=\frac{\hat{c}}{R_{p}}.


By requiring that :math:`\Delta t=\hat{\omega}_{p,0}^{-1}`, the electron
charge and mass can be fixed and calculated via

.. math::
|\hat{q}|=\frac{\hat{\omega}_{p,0}^{2}<\gamma>}{N_{ppc}\left(1+\frac{m_{-}}{m_{+}}\right)}.

The effect of this selection is that a plasma with the fiducial density :math:`n_0`, is always resolved numerically.
The skin depth is resolved by :math:`R_p` cells and the plasma oscillation time with :math:`1/\hat{c}` time steps.

Moving on to the relativistic cyclotron frequency for a given species

.. math::
\omega_{B}=\frac{qB}{mc\gamma},


to which the corresponding expression in code units can be found by
substitution with the appropriate scalings, e.g. :math:`B=B_{0}\hat{B}`.
This gives

.. math::
\omega_{B}=\frac{\hat{q}q_{0}}{\hat{m}m_{0}c\gamma}\frac{m_{0}c}{q_{0}\hat{c}\Delta t}\hat{B},


which cancels out and simplifies to

.. math::
\omega_{B}\Delta t=\frac{\hat{q}\hat{B}}{\hat{m}\hat{c}\gamma}\frac{m_{-}}{m_{+}},


where :math:`\hat{m}\frac{m_{-}}{m_{+}}`is the numerical particle mass
`1\hat{m}`for electron/positron plasma and :math:`1836\hat{m}` for protons.
Substituting for :math:`\Delta t` yields an expression for code unit cyclotron
frequency at any given scale :math:`\Delta x` and :math:`c`

.. math::
\omega_{B}\frac{\Delta x\hat{c}}{c}=\frac{\hat{q}\hat{B}}{\hat{m}\hat{c}\gamma}\frac{m_{-}}{m_{+}},


For the standard case where grid spacing and light speed are set to
one, the cyclotron frequency for the simulation can be calculated
via

.. math::
\omega_{B}=\frac{\hat{q}\hat{B}}{\hat{m}\hat{c}^{2}\gamma}\frac{m_{-}}{m_{+}},


Finally, the total plasma magnetization in Gaussian vs. code units
become

.. math::
\sigma=\frac{B^{2}}{4\pi nm_{e}c{{}^2}<\gamma>}=\left(\frac{\omega_{B}}{\omega_{p}}\right)^{2}=\frac{\hat{B}^{2}}{N_{ppc}\hat{m}\hat{c}^{2}<\gamma>}\left(1+\frac{m_{-}}{m_{+}}\right)^{-1}.

A given set of equations can be converted into the runko unit system using the scalings and definitions outlined above. In general, the choice of scalings, length and time scale will have terms cancelling out neatly and an expression will change from Gaussian units in the following ways when converting:

    i) Remove :math:`4\pi` Gaussian unit EM factors, 
    ii) Substitute physical variables with code versions, i.e. :math:`x\rightarrow \hat{x}` and so on,
    iii) Cancel out factors of :math:`q_e/m_e` since in this unit system :math:`|\hat{q}| = \hat{m}`.


Leapfrog Boris in Runko
=======================

The form of the relativistic Boris integrator must be changed to match
the unit system outlined above. Starting from the differenced version
of the Lorentz acceleration in code units (\ref{eq:lorentz_update})
for a particle with four-velocity :math:`\mathbf{u}(t^{n})=\mathbf{u}^{n}`

.. math::
\Delta[\mathbf{u}^{n}]=\frac{\hat{q}}{\hat{m}}\left(\hat{\mathbf{E}}+\frac{\hat{\mathbf{v}}}{\hat{c}}\times\hat{\mathbf{B}}\right).


The four-velocity :math:`\mathbf{u}^{n}` of is defined in units of :math:`c`
in runko to facilitate typical desired analysis, which require velocity
in these units. As the equations above define velocity in units of
`\hat{c}`, the velocity as tracked by runko must be translated to
and from these units for all calculations involving velocity. Velocities
with units of :math:`\hat{c}` rather than :math:`c` will be marked with primes,
and the conversion is simply a multiplication by :math:`\hat{c}`, so that
we have :math:`\mathbf{u}_{0}^{'}=\hat{c}\mathbf{u}^{n}`.

The first half of the electric field acceleration then take a form
similar to the standard formulation of Boris

.. math::
\mathbf{u}_{1}^{'}=\hat{c}\mathbf{u}_{0}+\frac{1}{2}\hat{s}\hat{\mathbf{E}},


but note the lack of an explicit time-step which has been cancelled
out by the choice of field scaling :math:`B_{0}`.
Here

.. math::
\hat{s}=\frac{\hat{q}}{\hat{m}_{s}m_{\pm}}=\frac{\mathrm{sign}\{\hat{q}_{0}\}}{m_{\pm}},


where :math:`m_{\pm}=m/m_{e}` is the particle mass in units of the electron rest mass. The Lorentz factor corresponding to the four-velocity :math:`\mathbf{u}_{1}^{'}` with units of :math:`\hat{c}` is given by

.. math::
\gamma_{1}=\frac{\sqrt{\hat{c}^{2}+\mathbf{u}_{1}^{'2}}}{\hat{c}}=\sqrt{1+\mathbf{u}_{1}^{'2}/\hat{c}}.


Proceeding, the first magnetic half-rotation takes the form

.. math::
\mathbf{u}_{2}=\mathbf{u}_{1}^{'}f+\hat{s}\frac{\mathbf{u}_{1}^{'}}{\hat{c}\gamma}\times\frac{1}{2}\hat{\mathbf{B}}f,


with 

.. math::
f=\frac{2}{1+\left(\frac{\hat{\mathbf{B}}}{2\hat{c}\gamma_{1}}\right)^{2}}.


For the final step, the last half rotation and acceleration from both
magnetic and electric fields can be combined to determine the final
velocity

.. math::
\mathbf{u'}_{3}=\mathbf{u'}_{1}+\hat{s}\frac{\mathbf{u'}_{2}}{\hat{c}\gamma}\times\frac{1}{2}\hat{\mathbf{B}}+\frac{1}{2}\hat{s}\hat{\mathbf{E}},


translating to units of :math:`c` then gives the updated four-velocity
at time :math:`t^{n+1}`

.. math::
\mathbf{u}^{n+1}=\frac{\mathbf{u'}_{3}}{\hat{c}}.


With velocity known, the position update takes the form

.. math::
\mathbf{x}^{n+1}=\mathbf{x}^{n}+\frac{\hat{c}\mathbf{u}^{n}}{\gamma(\hat{c}\mathbf{u}^{n})},


where position and velocity are staggered in time so that :math:`\mathbf{u}^{n+1}`
coincides with :math:`\mathbf{x}^{n+1/2}`. Note again the translation of
:math:`\mathbf{u}^{n}` from units of :math:`c` to units of :math:`\hat{c}`.

Relativistic Velocity-Verlet
============================

As an alternative to leapfrog integration, Velocity-Verlet can be
defined for relativistic particle motion. In velocity-Verlet, both
velocity and position are defined at the integer time-steps. Concurrent
velocity/position can simplify simulation setup as well as data handling,
since the half-step staggering of velocity need not be accounted for.
This is particularly important when comparing simulations of varying
time-step.

Starting again from relativistic Newton-Lorentz system

.. math::
\frac{d\mathbf{x}}{dt}=\mathbf{v}=\mathbf{\frac{u}{\gamma}},


.. math::
m\frac{d\mathbf{u}}{dt}=m\frac{d\gamma\mathbf{v}}{dt}=\mathbf{F}(\mathbf{x},\mathbf{v}),


where 

.. math::
\mathbf{F}(\mathbf{x},\mathbf{v})=q\left(\mathbf{E}(\mathbf{x})+\frac{\mathbf{v}}{c}\times\mathbf{B}(\mathbf{x})\right)=q\left(\mathbf{E}(\mathbf{x})+\frac{\mathbf{u}}{c\gamma}\times\mathbf{B}(\mathbf{x})\right),


and the relativistic factor is given by

.. math::
\gamma=\sqrt{\frac{1}{1-(\mathbf{v}/c)^{2}}}=\sqrt{1+\left(\frac{\mathbf{u}}{c}\right)^{2}}.


To achieve second order convergence accuracy like in Leapfrog, centre-difference
estimates of velocity and acceleration are applied and contracted
to one time-step length giving

.. math::
\frac{\mathbf{u}^{n+1/2}}{\gamma^{n+1/2}}=\frac{\mathbf{x}^{n+1}-\mathbf{x}^{n}}{\Delta t},


.. math::
\frac{1}{m}\mathbf{F}(\mathbf{x}^{n+1/2},\mathbf{v}^{n+1/2})=\frac{\mathbf{u}^{n+1}-\mathbf{u}^{n}}{\Delta t},


Rearranging these into update form for position and velocity yields

.. math::
\mathbf{x}^{n+1}=\mathbf{x}^{n}+\Delta t\frac{\mathbf{u}^{n+1/2}}{\gamma^{n+1/2}},


.. math::
\mathbf{u}^{n+1}=\mathbf{u}^{n}+\Delta t\frac{1}{m}\mathbf{F}(\mathbf{x}^{n+1/2},\frac{\mathbf{u}^{n+1/2}}{\gamma^{n+1/2}}).


The half-step velocity at :math:`t^{n+1/2}` in the position update can
be estimated by half an explicit Euler step of (\ref{eq:lorentz_rel})

.. math::
\frac{\mathbf{u}^{n+1/2}-\mathbf{u}^{n}}{\frac{1}{2}\Delta t}=\frac{1}{m}\mathbf{F}(\mathbf{x}^{n},\frac{\mathbf{u}^{n}}{\gamma^{n}}),


.. math::
\mathbf{u}^{n+1/2}=\mathbf{u}^{n}+\frac{\Delta t}{2}\frac{1}{m}\mathbf{F}(\mathbf{x}^{n},\frac{\mathbf{u}^{n}}{\gamma^{n}}),


The final system to be solved in cgs units is then is then

.. math::
\mathbf{u}^{n+1/2}=\mathbf{u}^{n}+\frac{\Delta t}{2}\mathbf{F}(\mathbf{x}^{n},\frac{\mathbf{u}^{n}}{\gamma^{n}}),


.. math::
\mathbf{x}^{n+1}=\mathbf{x}^{n}+\Delta t\frac{\mathbf{u}^{n+1/2}}{\gamma^{n+1/2}},


.. math::
\mathbf{u}^{n+1}=\mathbf{u}^{n}+\frac{\Delta t}{2}\left(\mathbf{F}(\mathbf{x}^{n},\frac{\mathbf{u}^{n}}{\gamma^{n}})+\mathbf{F}(\mathbf{x}^{n+1},\frac{\mathbf{u}^{n+1}}{\gamma^{n+1}})\right),


where all relativistic factors can be calculated using the corresponding
four-velocity and (). Note how the mass term has been absorbed by
`\mathbf{F}` to form the acceleration function

.. math::
\mathbf{F}(\mathbf{x},\frac{\mathbf{u}}{\gamma})=\frac{q}{m}\left(\mathbf{E}(\mathbf{x})+\frac{\mathbf{u}}{c\gamma}\times\mathbf{B}(\mathbf{x})\right).


The velocity update (\ref{eq:relvv_u}) is of the same form as in
leapfrog and can be solved explicitly using the Boris algorithm.

Relativistic Velocity-Verlet in runko
=====================================

In the section on leapfrog integration in runko above, the Boris algorithm
was derived for the specific case where runko scaling is in use. For
simple update equations however, there is no need to re-derive the
equations from the ground up using the correctly scaled variables.
As a worked example, the relativistic velocity-Verlet integrator will
be scaled for runko implementation below. Note that the unit-change
from :math:`c` to :math:`\hat{c}` of velocity will be omitted, as these are
specific to the storage scheme rather than scaling. To implement the
equations derived below, the velocities must simply be multiplied/divided
by :math:`\hat{c}` when retrieved and written to memory respectively.

To provide further context on the consequences of the normalisation
scheme used, it will also be demonstrated how the effective time-step
can be varied despite time-step being used to normalise the scheme. 

Substituting runko scalings into the relativistic velocity-verlet
equations (\ref{eq:relvv_uhalf})-(\ref{eq:relvv_u}) gives

.. math::
\mathbf{u}^{n+1/2}\frac{\text{\ensuremath{\Delta x}}}{\Delta t}=\mathbf{u}^{n}\frac{\text{\ensuremath{\Delta x}}}{\Delta t}+\frac{k_{\Delta t}\Delta t}{2}\mathbf{F}(\mathbf{\hat{x}}^{n}\Delta x,\frac{\mathbf{\hat{u}}^{n}}{\gamma^{n}}\frac{\Delta x}{\Delta t}),


.. math::
\mathbf{\hat{x}}^{n+1}\text{\ensuremath{\Delta x}}=\mathbf{\hat{x}}^{n}\text{\ensuremath{\Delta x}}+k_{\Delta t}\Delta t\frac{\mathbf{\hat{u}}^{n+1/2}}{\gamma^{n+1/2}}\frac{\text{\ensuremath{\Delta x}}}{\Delta t},


.. math::
\mathbf{\hat{u}}^{n+1}\frac{\Delta x}{\Delta t}=\mathbf{\hat{u}}^{n}\frac{\Delta x}{\Delta t}+\frac{k_{\Delta t}\Delta t}{2}\left(\mathbf{F}(\mathbf{\hat{x}}^{n}\Delta x,\frac{\mathbf{\hat{u}}^{n}}{\gamma^{n}}\frac{\Delta x}{\Delta t})+\mathbf{F}(\mathbf{\hat{x}}^{n+1}\Delta x,\frac{\mathbf{\hat{u}}^{n+1}}{\gamma^{n+1}}\frac{\Delta x}{\Delta t})\right).


The trick here is to recognise the difference between :math:`\Delta t`
as used to set the size of the time-step versus where it is used to
normalise the scheme. In the former case, time-step is rewritten as
a non-dimensional scale :math:`k_{\Delta t}` multiplied by the scaling
time-step :math:`\Delta t`. Carrying the factor :math:`k_{\Delta t}` through
the normalisation process allows simulations of varying effective
time-step to be performed for the same physical setup. The effective
time-step :math:`k_{\Delta t}` becomes a factor of a common normalisation
time-step :math:`\Delta t` between simulations. When :math:`k_{\Delta t}=1`
the simulation corresponds to standard runko \red{notation}. The same physics can
now be studied for varying temporal resolutions, by changing the number
of time-steps and ensuring a reciprocal change in effective time-step
`k_{\Delta t}`. For instance, when the number of time-steps is doubled,
the effective time-step must be halved.

Continuing the scaling process, expanding out :math:`\mathbf{F}`in velocity
half-step gives

.. math::
\mathbf{u}^{n+1/2}\frac{\text{\ensuremath{\Delta x}}}{\Delta t}=\mathbf{u}^{n}\frac{\text{\ensuremath{\Delta x}}}{\Delta t}+\frac{k_{\Delta t}\Delta t}{2}\frac{q_{0}\hat{q}}{m_{0}\hat{m}}B_{0}\left(\mathbf{\hat{E}}^{n}+\frac{\mathbf{\hat{u}}^{n}}{c\gamma^{n}}\frac{\Delta x}{\Delta t}\times\hat{\mathbf{B}}^{n}\right),


.. math::
\mathbf{u}^{n+1/2}=\mathbf{u}^{n}+\frac{\text{\ensuremath{\Delta t}}}{\Delta x}\frac{k_{\Delta t}\Delta t}{2}\frac{q_{0}\hat{q}}{m_{0}\hat{m}}\frac{m_{0}c}{q_{0}\hat{c}\Delta t}\left(\mathbf{\hat{E}}^{n}+\frac{\mathbf{\hat{u}}^{n}}{c\gamma^{n}}\frac{\Delta x}{\Delta t}\times\hat{\mathbf{B}}^{n}\right),


Remembering the relation :math:`\Delta t/\Delta x=\hat{c}/c` and cancelling
out terms accordingly in both the position update and velocity half-step,
these equations must take the following form in runko

.. math::
\mathbf{u}^{n+1/2}=\mathbf{u}^{n}+\frac{k_{\Delta t}}{2}\frac{\hat{q}}{\hat{m}}\left(\mathbf{\hat{E}}^{n}+\frac{\mathbf{\hat{u}}^{n}}{\hat{c}\gamma^{n}}\times\hat{\mathbf{B}}^{n}\right),


.. math::
\mathbf{\hat{x}}^{n+1}=\mathbf{\hat{x}}^{n}+k_{\Delta t}\frac{\mathbf{\hat{u}}^{n+1/2}}{\gamma^{n+1/2}},


Meanwhile, expanding out :math:`\mathbf{F}`in velocity update yields

.. math::
\mathbf{\hat{u}}^{n+1}\frac{\Delta x}{\Delta t}=\mathbf{\hat{u}}^{n}\frac{\Delta x}{\Delta t}+\frac{k_{\Delta t}\Delta t}{2}\left(\frac{q_{0}\hat{q}}{m_{0}\hat{m}}B_{0}\left(\mathbf{\hat{E}}^{n}+\frac{\hat{\mathbf{v}}^{n}}{c}\frac{c}{\hat{c}}\times\hat{\mathbf{B}}^{n}\right)+\frac{q_{0}\hat{q}}{m_{0}\hat{m}}B_{0}\left(\mathbf{\hat{E}}^{n+1}+\frac{\hat{\mathbf{v}}^{n+1}}{c}\frac{c}{\hat{c}}\times\hat{\mathbf{B}}^{n+1}\right)\right),


.. math::
\mathbf{\hat{u}}^{n+1}\frac{\Delta x}{\Delta t}=\mathbf{\hat{u}}^{n}\frac{\Delta x}{\Delta t}+\frac{k_{\Delta t}\Delta t}{2}\frac{q_{0}\hat{q}}{m_{0}\hat{m}}B_{0}\left(\left(\mathbf{\hat{E}}^{n}+\frac{\hat{\mathbf{v}}^{n}}{\hat{c}}\times\hat{\mathbf{B}}^{n}\right)+\left(\mathbf{\hat{E}}^{n+1}+\frac{\hat{\mathbf{v}}^{n+1}}{\hat{c}}\times\hat{\mathbf{B}}^{n+1}\right)\right),


.. math::
\mathbf{\hat{u}}^{n+1}\frac{\Delta x}{\Delta t}=\mathbf{\hat{u}}^{n}\frac{\Delta x}{\Delta t}+\frac{k_{\Delta t}\Delta t}{2}\frac{q_{0}\hat{q}}{m_{0}\hat{m}}\frac{m_{0}c}{q_{0}\hat{c}\Delta t}\left(\mathbf{\hat{E}}^{n}+\mathbf{\hat{E}}^{n+1}+\frac{\hat{\mathbf{v}}^{n}}{\hat{c}}\times\hat{\mathbf{B}}^{n}+\frac{\hat{\mathbf{v}}^{n+1}}{\hat{c}}\times\hat{\mathbf{B}}^{n+1}\right),


.. math::
\mathbf{\hat{u}}^{n+1}=\mathbf{\hat{u}}^{n}+k_{\Delta t}\frac{\hat{q}}{\hat{m}}\left(\frac{\mathbf{\hat{E}}^{n}+\mathbf{\hat{E}}^{n+1}}{2}+\frac{\hat{\mathbf{v}}^{n}+\hat{\mathbf{v}}^{n+1}}{2\hat{c}}\times\frac{\hat{\mathbf{B}}^{n}+\hat{\mathbf{B}}^{n+1}}{2}\right),


where :math:`\hat{\mathbf{v}}^{n}=\hat{\mathbf{u}}^{n}/\gamma^{n}`.

The seemingly implicit form of the velocity update (\ref{eq:relvv_u_runko})
is of the same form as the leapfrog discretisation of velocity and
can now be solved with the runko-scaled Boris algorithm above. Together
with (\ref{eq:relvv_uhalf_runko}) and (\ref{eq:relvv_x_runko}) the
evolution in time of particle motion can now be solved to second order,
so long as :math:`\mathbf{\hat{E}}(\mathbf{\hat{x}},\mathbf{\hat{u}})`
and :math:`\mathbf{\hat{B}}(\mathbf{\hat{x}},\mathbf{\hat{u}})` can be
evaluated.



