# Equations of Motion: Mathematical Derivation

## 1. Newton's Law of Universal Gravitation

The gravitational force on body $i$ due to body $j$ is:

$$\mathbf{F}_{ij} = G \frac{m_i m_j}{|\mathbf{r}_j - \mathbf{r}_i|^2} \hat{\mathbf{r}}_{ij} = G \frac{m_i m_j (\mathbf{r}_j - \mathbf{r}_i)}{|\mathbf{r}_j - \mathbf{r}_i|^3}$$

**Newton's Third Law:** $\mathbf{F}_{ij} = -\mathbf{F}_{ji}$ — forces are equal and opposite.

## 2. Three-Body Equations of Motion

Applying Newton's second law $\mathbf{F} = m\mathbf{a}$ with all pairwise forces:

$$m_1 \ddot{\mathbf{r}}_1 = G\frac{m_1 m_2(\mathbf{r}_2 - \mathbf{r}_1)}{|\mathbf{r}_{12}|^3} + G\frac{m_1 m_3(\mathbf{r}_3 - \mathbf{r}_1)}{|\mathbf{r}_{13}|^3}$$

$$m_2 \ddot{\mathbf{r}}_2 = G\frac{m_1 m_2(\mathbf{r}_1 - \mathbf{r}_2)}{|\mathbf{r}_{12}|^3} + G\frac{m_2 m_3(\mathbf{r}_3 - \mathbf{r}_2)}{|\mathbf{r}_{23}|^3}$$

$$m_3 \ddot{\mathbf{r}}_3 = G\frac{m_1 m_3(\mathbf{r}_1 - \mathbf{r}_3)}{|\mathbf{r}_{13}|^3} + G\frac{m_2 m_3(\mathbf{r}_2 - \mathbf{r}_3)}{|\mathbf{r}_{23}|^3}$$

Dividing each by $m_i$ to get accelerations:

$$\ddot{\mathbf{r}}_i = G \sum_{j \neq i} \frac{m_j (\mathbf{r}_j - \mathbf{r}_i)}{|\mathbf{r}_j - \mathbf{r}_i|^3}$$

## 3. Hamiltonian Formulation

The system has a Hamiltonian $H = T + V$:

$$T = \frac{1}{2}\sum_i m_i |\dot{\mathbf{r}}_i|^2 \qquad (\text{kinetic energy})$$

$$V = -G\left(\frac{m_1 m_2}{r_{12}} + \frac{m_1 m_3}{r_{13}} + \frac{m_2 m_3}{r_{23}}\right) \qquad (\text{potential energy})$$

Hamilton's equations:

$$\dot{\mathbf{r}}_i = \frac{\partial H}{\partial \mathbf{p}_i} = \frac{\mathbf{p}_i}{m_i}, \qquad \dot{\mathbf{p}}_i = -\frac{\partial H}{\partial \mathbf{r}_i} = \mathbf{F}_i$$

## 4. Conservation Laws from Symmetry (Noether's Theorem)

| Symmetry | Conserved Quantity | Expression |
|----------|------------------|-----------|
| Time translation ($t \to t + \varepsilon$) | Total energy | $E = T + V$ |
| Spatial translation ($\mathbf{r} \to \mathbf{r} + \boldsymbol{\varepsilon}$) | Linear momentum | $\mathbf{P} = \sum_i m_i \dot{\mathbf{r}}_i$ |
| Rotation ($\theta \to \theta + \varepsilon$) | Angular momentum | $\mathbf{L} = \sum_i m_i (\mathbf{r}_i \times \dot{\mathbf{r}}_i)$ |

## 5. Center of Mass Frame

Define the center of mass:
$$\mathbf{R}_\text{cm} = \frac{\sum_i m_i \mathbf{r}_i}{\sum_i m_i} = \frac{m_1\mathbf{r}_1 + m_2\mathbf{r}_2 + m_3\mathbf{r}_3}{M}$$

Since $\sum_i \mathbf{F}_i = 0$ (Newton's third law), $\ddot{\mathbf{R}}_\text{cm} = 0$, meaning the CM moves at constant velocity.

**In the CM frame:** $\mathbf{R}_\text{cm} = \mathbf{0}$ and $\dot{\mathbf{R}}_\text{cm} = \mathbf{0}$ at all times.

Transforming to the CM frame eliminates 4 of the 12 degrees of freedom.

## 6. Softening Regularization

The gravitational force diverges as $r_{ij} \to 0$. To regularize:

$$|\mathbf{r}_{ij}|^3 \to \left(|\mathbf{r}_{ij}|^2 + \varepsilon^2\right)^{3/2}$$

For $r \gg \varepsilon$: softened force ≈ exact force
For $r \ll \varepsilon$: force remains finite instead of diverging

Default: $\varepsilon = 10^{-4}$ AU.

## 7. Units and Dimensional Analysis

Using **astronomical units** where $G = 4\pi^2$:

$$G = 4\pi^2 \;\frac{\text{AU}^3}{M_\odot \cdot \text{yr}^2}$$

**Verification:** Earth orbits Sun at 1 AU with period 1 yr. Kepler's third law:
$$T^2 = \frac{4\pi^2}{GM} a^3 \;\Rightarrow\; 1 = \frac{4\pi^2}{4\pi^2 \cdot 1} \cdot 1 = 1 \quad \checkmark$$

## 8. Degrees of Freedom Count

In 2D:
- 3 bodies × 2 position components = 6 position DOF
- 3 bodies × 2 velocity components = 6 velocity DOF
- Total: **12-dimensional phase space**

Conserved quantities reduce effective DOF:
- Energy: −1 DOF
- Angular momentum: −1 DOF
- Linear momentum (CM frame): −2 DOF

**Effective phase space dimension: 8** — still high-dimensional enough for chaos.

## 9. Connection to General Relativity

The Newtonian geodesic equation:

$$\ddot{\mathbf{r}}_i = -\nabla \Phi(\mathbf{r}_i) \qquad \text{where } \Phi = -G\sum_{j\neq i}\frac{m_j}{|\mathbf{r}_j - \mathbf{r}_i|}$$

The GR geodesic equation in curved spacetime:

$$\frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\alpha\beta} \frac{dx^\alpha}{d\lambda}\frac{dx^\beta}{d\lambda} = 0$$

The Christoffel symbols $\Gamma^\mu_{\alpha\beta}$ play the role of $\nabla\Phi$ — they encode how spacetime curvature (gravity) accelerates free particles. In the weak-field, slow-motion limit, the GR geodesic equation reduces exactly to the Newtonian equation above.