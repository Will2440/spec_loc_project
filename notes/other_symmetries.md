

### 1. Local Orbital / Sublattice Symmetry ($\sigma_z$)

* **The Concept:** Instead of acting on the spatial coordinates of the lattice, this operator acts purely on the internal degrees of freedom (like the $s$ and $p$ orbitals, or $A$ and $B$ sublattices).
* 
**The Motivation:** Global spatial symmetries violently anti-commute with the position operators $X$ and $Y$. A purely internal, local symmetry perfectly commutes with the real-space coordinates.



### 2. Discrete Translation / Momentum Filtering ($T$)

* **The Concept:** Utilizing the underlying periodicity of the lattice itself.
* 
**The Motivation:** This is the only spatial symmetry that strictly separates the $\Gamma(0,0)$ point from the $M(\pi,\pi)$ point mathematically. A state at $\Gamma$ has a uniform real-space phase, while a state at $M$ acts as a spatial checkerboard. We could attempt to use a translation operator or a spatial envelope filter to target these phase gradients.



### 3. Point Group Symmetries (Rotations & Mirrors)

* **The Concept:** Exploiting the specific geometric symmetries of the square lattice, such as four-fold rotation ($C_4$) or mirror reflection planes.
* 
**The Motivation:** By looking at the "little group" of symmetries that leave specific $\mathbf{k}$ points invariant (or map them to one another), we might find a geometric operator that distinguishes the $\Gamma$ pocket from the $M$ pocket in a way that inversion could not.



### 4. Anti-Unitary Symmetries (Time-Reversal & Particle-Hole)

* **The Concept:** The symmetries that define the Altland-Zirnbauer topological classification table.
* 
**The Motivation:** While we know the standard QWZ model explicitly breaks Time-Reversal Symmetry to achieve a Chern insulator phase , and generically violates Particle-Hole symmetry, we must formally evaluate if defining a chiral proxy or leveraging a pseudo-symmetry could aid in state separation.



---

To evaluate these rigorously, we need to subject each one to the same three tests: Does it separate the states in momentum space? Does it commute with the real-space localizer? Does it survive physical disorder?