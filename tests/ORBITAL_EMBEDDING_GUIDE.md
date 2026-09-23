# Orbital Sublattice Embedding Implementation

## Overview
Implementation of k-space orbital position tuning via phase shift transformation in the QWZ Hamiltonian. This allows controlled adjustment of orbital sublattice embedding without changing the base model parameters.

## Mathematical Foundation

### Phase Shift Transformation
The orbital embedding is controlled by a unitary transformation of the Hamiltonian:

$$H'(\mathbf{k}) = e^{-i\theta(\mathbf{k},\phi) \sigma_z} H(\mathbf{k}) e^{i\theta(\mathbf{k},\phi) \sigma_z}$$

where the phase angle is:
$$\theta(\mathbf{k},\phi) = \mathbf{k} \cdot \mathbf{d}(\phi) = k_x d_x(\phi) + k_y d_y(\phi)$$

### Displacement Vector
The orbital displacement vector is parameterized by angle `phi` and magnitude `orbital_displacement`:

$$\mathbf{d}(\phi) = d_0 \begin{pmatrix} \cos(\phi) \\ \sin(\phi) \end{pmatrix}$$

where:
- `d₀ = orbital_displacement`: Controls the magnitude of displacement
- `φ`: Controls the direction (angular position) in k-space

### Hamiltonian Transformation
Under the $\sigma_z$-based rotation, the Pauli matrices transform as:

$$\sigma_x \to \cos(2\theta) \sigma_x + \sin(2\theta) \sigma_y$$
$$\sigma_y \to -\sin(2\theta) \sigma_x + \cos(2\theta) \sigma_y$$
$$\sigma_z \to \sigma_z \text{ (invariant)}$$

This gives the rotated d-vector components:
$$d'_x = d_x \cos(2\theta) - d_y \sin(2\theta)$$
$$d'_y = d_x \sin(2\theta) + d_y \cos(2\theta)$$

## API Changes

### `k_space_perturbed_qwz_hamiltonian`
**New parameters:**
- `phi::Real = 0.0` — Angular position of displacement vector (radians)
- `orbital_displacement::Real = 0.0` — Magnitude of displacement |d|

**Example:**
```julia
H = k_space_perturbed_qwz_hamiltonian(
    kx, ky;
    A=1.0, B=1.0, m=-1.0,
    phi=π/4,            # 45° direction
    orbital_displacement=0.2  # Magnitude
)
```

### `compute_bulk_band_berry_data`
Now passes through `phi` and `orbital_displacement` parameters:

```julia
energies = compute_bulk_band_berry_data(
    Nkx=101, Nky=101,
    A=1.0, B=1.0, m=-1.0,
    phi=π/6,
    orbital_displacement=0.15
)
```

## Physical Interpretation

### Effect on Band Structure
- `phi = 0`: Displacement along kₓ-direction → modifies x-hopping
- `phi = π/2`: Displacement along k_y-direction → modifies y-hopping  
- `phi = π/4`: Diagonal displacement → mixes both directions
- `orbital_displacement`: Strength of the embedding effect (0 = no effect)

### Use Cases
1. **Tuning orbital hybridization**: Adjust sublattice orbital overlap
2. **Controlling band gap**: Varies with both `phi` and `orbital_displacement`
3. **Engineering topology**: Alters Berry curvature and Chern numbers
4. **Simulating disorder**: Can represent off-site perturbations

## Implementation Details

The transformation is applied only when `|orbital_displacement| > 1e-14` to avoid numerical issues with near-zero values.

**Phase shift angle calculation:**
```julia
theta = kx * (orbital_displacement * cos(phi)) + 
        ky * (orbital_displacement * sin(phi))
```

**Rotation application:**
```julia
cos_2theta = cos(2 * theta)
sin_2theta = sin(2 * theta)

d_x_rot = d_x * cos_2theta - d_y * sin_2theta
d_y_rot = d_x * sin_2theta + d_y * cos_2theta
```

## Example Scan

To scan the effect of orbital embedding:

```julia
phi_vals = range(0, 2π; length=16)
d_vals = range(0, 0.5; length=11)

for phi in phi_vals
    for d in d_vals
        H = k_space_perturbed_qwz_hamiltonian(
            kx, ky;
            phi=phi, 
            orbital_displacement=d
        )
        # Compute properties...
    end
end
```

## Combination with Other Perturbations

The phase shift transformation is applied **after** winding distortions but **before** the final Hamiltonian assembly, allowing simultaneous use with:
- Symmetric perturbations (`perturbation_type=:symmetric`)
- Tilt perturbations (`perturbation_type=:tilt`)
- Winding distortions (`winding_number > 0`)

## Reversibility

The transformation is unitary and therefore reversible:
- Setting `orbital_displacement=0` recovers the original QWZ Hamiltonian
- The transformation preserves eigenvalue spectrum structure
