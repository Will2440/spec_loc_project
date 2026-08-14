# Pedagogical Summary: A Band-Graded Spectral Localiser for Indirect-Gap Topology

## 1. The problem with the ordinary spectral localiser

For a two-dimensional Hamiltonian (H), the ordinary spectral localiser at position ((x_0,y_0)), energy (E), and spatial-energy scaling parameter (\kappa) is

[
L_\kappa(H;E,x_0,y_0)
=====================

\begin{pmatrix}
H-E &
\kappa\bigl((X-x_0)-i(Y-y_0)\bigr)
[2mm]
\kappa\bigl((X-x_0)+i(Y-y_0)\bigr) &
-(H-E)
\end{pmatrix}.
]

Equivalently,

[
L_\kappa
========

\kappa(X-x_0)\otimes \sigma_x
+
\kappa(Y-y_0)\otimes \sigma_y
+
(H-E)\otimes \sigma_z .
]

When (E) lies in a spectral gap of (H), the half-signature

[
\operatorname{Ind}_{\mathrm{loc}}(H;E)
======================================

\frac{1}{2}\operatorname{sig}L_\kappa
]

can reproduce the Chern number, provided the finite system is sufficiently large and (\kappa) lies in an appropriate interval.

The key hypothesis is the existence of an energy gap:

[
(H-E)^2\geq g_E^2 I
]

for some (g_E>0). Squaring the localiser gives, schematically,

[
L_\kappa^2
==========

(H-E)^2
+
\kappa^2\bigl((X-x_0)^2+(Y-y_0)^2\bigr)
+
\text{commutator terms}.
]

The positive contribution ((H-E)^2\geq g_E^2I) must dominate the commutator errors. A typical sufficient condition is therefore of the form

[
2\kappa M_H<g_E^2,
\qquad
M_H=
\max{|[X,H]|,|[Y,H]|}.
]

This explains both the strength and the weakness of the ordinary localiser:

* it works directly with the physical Hamiltonian;
* it resolves topology relative to a chosen energy;
* but it needs a spectral or local spectral gap at that energy.

In the indirectly overlapping QWZ model,

[
H_\gamma(\mathbf k)
===================

H_{\mathrm{QWZ}}(\mathbf k)
+
\gamma(\cos k_x+\cos k_y)I_2,
]

the identity-valued perturbation changes the two band energies without changing their eigenvectors or eigenprojectors. Thus the direct separation between the two bands can remain open at every (\mathbf k), while the energy ranges of the two bands overlap globally. Once this happens, no energy (E) separates all lower-band states from all upper-band states. The ordinary localiser consequently loses the energy gap that supported its invertibility, even though the lower-band bundle and its Chern number remain unchanged.

The failure is therefore not that the topology has disappeared. The failure is that **energy ordering no longer identifies the band carrying the topology**.

---

## 2. Why ordinary spectral flattening cannot repair the problem

In a spectrally gapped insulator one often forms the flattened Hamiltonian

[
\operatorname{sgn}(H-E_F)
=========================

I-2\chi(H<E_F).
]

This is an energy flattening: every state below (E_F) is assigned (-1), while every state above (E_F) is assigned (+1).

That construction is unavailable after indirect band overlap. If a lower-band state and an upper-band state have the same energy, then every scalar function (f(H)) acts on them using the same value (f(E)). It therefore cannot assign one state to the lower band and the other to the upper band.

Consequently, there is no scalar function (f) for which

[
f(H_\gamma)=P_-
]

once the two band spectra overlap. In particular, the desired projector cannot be obtained from:

[
\Theta(E_F-H_\gamma),
\qquad
\operatorname{sgn}(H_\gamma-E_F),
]

a Fermi–Dirac filter, or a KPM approximation to any scalar energy filter. The proof document identifies this as the spectral-function no-go theorem.

The repaired construction must therefore use **band information in addition to the unordered spectrum of (H)**.

---

## 3. Replace energy flattening by band flattening

Let (P_-) denote the projector onto the chosen lower-band bundle and (P_+=I-P_-). Define the band grading

[
Q=P_+-P_-=I-2P_-.
]

This operator satisfies

[
Q^\dagger=Q,
\qquad
Q^2=I,
\qquad
\sigma(Q)={-1,+1}.
]

The negative spectral subspace of (Q) is exactly the desired lower-band subspace:

[
\chi(Q<0)=P_-.
]

The central distinction is

[
\boxed{
\text{energy flattening } \operatorname{sgn}(H-E_F)
\quad\neq\quad
\text{band flattening } Q=I-2P_-.
}
]

Energy flattening asks:

> Is the energy of this state above or below (E_F)?

Band flattening asks:

> Does this state belong to the selected lower or upper band?

These questions are equivalent in a globally gapped two-band insulator, but they cease to be equivalent after indirect overlap.

For the clean perturbed QWZ model, write

[
H_\gamma(\mathbf k)
===================

d_0(\mathbf k)I_2
+
\mathbf d(\mathbf k)\cdot\boldsymbol{\sigma},
]

where

[
d_0(\mathbf k)
==============

\gamma(\cos k_x+\cos k_y).
]

The band grading is

[
Q_0(\mathbf k)
==============

\widehat{\mathbf d}(\mathbf k)\cdot\boldsymbol{\sigma},
\qquad
\widehat{\mathbf d}
===================

\frac{\mathbf d}{|\mathbf d|}.
]

The scalar term (d_0I_2) does not enter (Q_0). Hence

[
Q_{0,\gamma}(\mathbf k)=Q_{0,0}(\mathbf k)
]

for every (\gamma). This is why the new construction is insensitive to the strength of the indirect overlap in the clean toy model.

---

## 4. The band-graded spectral localiser

The repaired localiser is obtained by replacing the physical energy slot (H-E) by the band grading (Q):

[
\widetilde L_\kappa(Q;x_0,y_0)
==============================

\begin{pmatrix}
Q &
\kappa\bigl((X-x_0)-i(Y-y_0)\bigr)
[2mm]
\kappa\bigl((X-x_0)+i(Y-y_0)\bigr) &
-Q
\end{pmatrix}.
]

Equivalently,

[
\widetilde L_\kappa
===================

\kappa(X-x_0)\otimes\sigma_x
+
\kappa(Y-y_0)\otimes\sigma_y
+
Q\otimes\sigma_z .
]

The proposed invariant is

[
C_{\mathrm{graded}}
===================

\frac{1}{2}\operatorname{sig}\widetilde L_\kappa .
]

The matrix structure is identical to that of the ordinary spectral localiser. The essential change is the operator placed in its third, or “energy”, slot:

[
\begin{array}{c|c|c}
&\text{ordinary localiser}&\text{band-graded localiser}\
\hline
\text{slot operator}
&H-E
&Q\
\text{reference value}
&E
&0\
\text{separates}
&\text{states by energy}
&\text{states by band}\
\text{gap used}
&\operatorname{dist}(E,\sigma(H))
&\operatorname{dist}(0,\sigma(Q))=1\
\text{negative projector}
&\chi(H<E)
&P_-\
\text{indirect overlap}
&\text{generally fails}
&\text{remains well defined if (Q) is available}
\end{array}
]

The spectral-localiser theorem does not fundamentally require the slot operator to be the physical Hamiltonian. It requires a bounded, self-adjoint, sufficiently local operator that is spectrally gapped at the reference value. The grading (Q) is precisely such an operator when it is local or quasi-local. Its negative projection is the band projector whose index is the Chern number.

---

## 5. Why the graded localiser retains a gap

Assume for the moment that (Q) is sufficiently local that

[
M_Q
===

\max{|[X,Q]|,|[Y,Q]|}
<\infty.
]

Squaring the graded localiser gives

[
\begin{aligned}
\widetilde L_\kappa^2
={}&
\left[
\kappa^2(X-x_0)^2+
\kappa^2(Y-y_0)^2+
Q^2
\right]\otimes I
\
&+
\kappa[X,Q]\otimes\sigma_x\sigma_z
+
\kappa[Y,Q]\otimes\sigma_y\sigma_z .
\end{aligned}
]

Because (Q^2=I),

[
\widetilde L_\kappa^2
\geq
\left(1-2\kappa M_Q\right)I.
]

Thus, whenever

[
\kappa<\frac{1}{2M_Q},
]

the localiser is invertible and has the lower gap estimate

[
\operatorname{gap}(\widetilde L_\kappa)
\geq
\sqrt{1-2\kappa M_Q}.
]

This is the basic mechanism of the repair.

For the ordinary localiser, the corresponding positive term is ((H-E)^2), so its lower bound is (g_E^2), which vanishes when no global spectral gap exists.

For the graded localiser, the positive term is

[
Q^2=I
]

exactly. The grading therefore supplies a fixed unit gap even when the spectrum of the physical Hamiltonian is strongly overlapping. The only competition is with the position commutators of (Q), which measure its nonlocality.

This also clarifies the role of locality. The method does not succeed merely because one has found an involution (Q). It succeeds because one has found a **local or quasi-local involution** whose negative subspace is the desired band.

---

## 6. Why the half-signature gives the Chern number

Let

[
P_-=\frac{I-Q}{2}.
]

For a local, gapped (Q), the even-dimensional spectral-localiser theorem gives

[
\frac{1}{2}\operatorname{sig}\widetilde L_{\kappa,\rho}
=======================================================

\operatorname{Ind}(P_-FP_-),
]

where (F) is the Dirac phase constructed from the position operators. The index on the right is the standard real-space index pairing for the projection (P_-), and equals its Chern number:

[
\operatorname{Ind}(P_-FP_-)=C(P_-).
]

The new method therefore does not define an unrelated invariant. It computes the usual Chern invariant of the selected band projector, but through a localiser whose gap is tied to the grading rather than to the physical energy spectrum.

In the clean perturbed QWZ model,

[
Q_\gamma=Q_0
]

for all (\gamma). Hence the following quantities should all be independent of (\gamma):

[
|[X,Q]|,
\qquad
|[Y,Q]|,
\qquad
\operatorname{sig}\widetilde L_\kappa,
\qquad
\operatorname{gap}(\widetilde L_\kappa).
]

Meanwhile, the ordinary localiser gap should deteriorate and eventually close as the physical bands begin to overlap indirectly.

This comparison is the most direct numerical illustration of what has changed:

> The ordinary localiser follows the energy dispersion; the graded localiser follows the band geometry.

---

## 7. The three computational regimes must be distinguished

### 7.1 Clean system

This is the fully justified test case.

The grading can be constructed directly in Bloch space:

[
Q_0(\mathbf k)
==============

\widehat{\mathbf d}(\mathbf k)\cdot\boldsymbol{\sigma}.
]

One then Fourier-transforms it to real space, restricts it to a finite box, constructs (\widetilde L_\kappa), and computes its half-signature.

The expected result is

[
\frac{1}{2}\operatorname{sig}\widetilde L_\kappa=C
]

for every (\gamma), including values for which the physical Hamiltonian has no global energy gap. This clean implementation is presently on firm theoretical ground.

### 7.2 Structured disorder

Suppose the Hamiltonian admits a decomposition

[
H=h+g
]

such that:

1. (h) has a spectral or mobility gap at zero;
2. (P_-=\chi(h<0));
3. the remaining component (g) commutes with the grading (Q=I-2P_-).

Then (Q) can be constructed from the gapped auxiliary operator (h), even when the complete physical Hamiltonian (H) is metallic.

For a genuine spectral gap of (h), one may calculate

[
P_-=\chi(h<0)
]

by diagonalisation or by a contour integral around the negative spectrum:

[
P_-=
\frac{1}{2\pi i}
\oint_{\mathcal C}(z-h)^{-1},dz.
]

This is not forbidden by the spectral-function no-go theorem, because the functional calculus is applied to the auxiliary gapped operator (h), not to the overlapping physical Hamiltonian (H).

### 7.3 Weak generic disorder

With generic disorder, the clean band grading (Q_0) will not normally commute with the perturbed Hamiltonian. One must instead search for an adapted grading (Q) near (Q_0).

The proved weak-disorder construction is based on a contraction or fixed-point argument. Numerically, this means starting from (Q_0), iteratively rotating away the part of (H) that is off-diagonal with respect to the current grading, and stopping when the residual is sufficiently small.

The important certificate is

[
|Q-Q_0|<1.
]

This guarantees that (Q) remains in the same homotopy class as the clean grading. The adapted grading, its quasi-locality, and its homotopy class are justified in the weak-disorder regime.

There remains an important qualification: the final finite-volume theorem identifying the disordered truncated localiser’s half-signature with the Chern number, including all volume-dependent constants, is still pending correction. Consequently, weak-disorder half-signatures should presently be labelled **numerically supported**, rather than cited as a completed finite-volume theorem.

Strong generic disorder deep in the overlap regime should not yet be treated as theoretically certified. The loss of fixed-point convergence is instead a useful numerical indicator of the boundary of the weak-disorder construction.

---

# Part II: Pseudo-Implementation

## 8. Common localiser routine

The ordinary and graded localisers should use exactly the same numerical routine. The routine should accept an arbitrary Hermitian slot operator (A):

[
A=
\begin{cases}
H-EI,&\text{ordinary localiser},\
Q,&\text{graded localiser}.
\end{cases}
]

This ensures that numerical differences arise from the change (H-E\mapsto Q), not from two different implementations.

```julia
"""
Construct the two-dimensional spectral localiser for a Hermitian
slot operator A.

Use:
    A = H - E*I    for the ordinary localiser
    A = Q          for the band-graded localiser
"""
function build_localiser(
    A::AbstractMatrix,
    X::AbstractMatrix,
    Y::AbstractMatrix,
    x0::Real,
    y0::Real;
    kappa::Real
)
    D = size(A, 1)

    @assert size(A) == (D, D)
    @assert size(X) == (D, D)
    @assert size(Y) == (D, D)

    I_D = Matrix{ComplexF64}(I, D, D)

    X0 = Matrix{ComplexF64}(X) - x0 * I_D
    Y0 = Matrix{ComplexF64}(Y) - y0 * I_D
    A0 = Matrix{ComplexF64}(A)

    L = [
        A0                         kappa * (X0 - im * Y0);
        kappa * (X0 + im * Y0)    -A0
    ]

    return Hermitian(L)
end
```

Then use one diagnostic routine for both cases:

```julia
"""
Return:
    signature = number of positive eigenvalues minus negative eigenvalues
    gap       = minimum absolute eigenvalue
    index     = signature / 2, if the localiser is numerically invertible
"""
function localiser_signature_and_gap(
    A::AbstractMatrix,
    X::AbstractMatrix,
    Y::AbstractMatrix,
    x0::Real,
    y0::Real;
    kappa::Real,
    zero_tol::Real = 1e-10
)
    L = build_localiser(A, X, Y, x0, y0; kappa=kappa)
    λ = eigvals(L)

    npos = count(x -> x >  zero_tol, λ)
    nneg = count(x -> x < -zero_tol, λ)

    signature = npos - nneg
    gap = minimum(abs, λ)

    index = gap > zero_tol ? signature / 2 : missing

    return (
        signature = signature,
        index = index,
        gap = gap
    )
end
```

For large systems, the full eigendecomposition should eventually be replaced by:

1. a Hermitian-indefinite factorisation to obtain the inertia and hence the signature;
2. a smallest-magnitude eigensolver to estimate the localiser gap.

The dense eigensolver should nevertheless be retained as the reference implementation for small systems. The numerical schema correctly recommends validating the faster inertia computation against full diagonalisation before relying on it.

---

## 9. Clean construction of (Q_0) in Bloch space

For a QWZ convention of the form

[
H_0(\mathbf k)
==============

d_0(\mathbf k)I+
d_x(\mathbf k)\sigma_x+
d_y(\mathbf k)\sigma_y+
d_z(\mathbf k)\sigma_z,
]

construct

[
Q_0(\mathbf k)
==============

\frac{
d_x(\mathbf k)\sigma_x+
d_y(\mathbf k)\sigma_y+
d_z(\mathbf k)\sigma_z
}{
\sqrt{d_x^2+d_y^2+d_z^2}
}.
]

The identity coefficient (d_0), including the (\gamma)-dependent indirect-gap perturbation, must be discarded.

```julia
function qwz_bloch_grading(
    kx::Real,
    ky::Real;
    A::Real,
    B::Real,
    m::Real
)
    # These components must match the real-space QWZ convention.
    dx = A * sin(kx)
    dy = A * sin(ky)
    dz = m + B * (2 - cos(kx) - cos(ky))

    dnorm = sqrt(dx^2 + dy^2 + dz^2)

    dnorm > 0 || error(
        "Direct band gap closes at k = ($kx, $ky)."
    )

    Qk = (
        dx * sigma_x +
        dy * sigma_y +
        dz * sigma_z
    ) / dnorm

    return ComplexF64.(Qk)
end
```

Notice that there is deliberately no `gamma` argument. A fundamental clean-system unit test is therefore

```julia
@assert Qk_gamma1 ≈ Qk_gamma2
```

for arbitrary values of (\gamma_1) and (\gamma_2).

---

## 10. Fourier transformation to a finite real-space grading

On a periodic (L_x\times L_y) momentum grid,

[
\mathbf k_{n_x,n_y}
===================

\left(
\frac{2\pi n_x}{L_x},
\frac{2\pi n_y}{L_y}
\right),
]

compute the real-space blocks

[
Q_{\mathbf r,\mathbf r'}
========================

\frac{1}{L_xL_y}
\sum_{\mathbf k}
e^{i\mathbf k\cdot(\mathbf r-\mathbf r')}
Q_0(\mathbf k).
]

A transparent reference implementation is:

```julia
function clean_grading_real_space(
    Lx::Int,
    Ly::Int;
    A::Real,
    B::Real,
    m::Real
)
    nsites = Lx * Ly
    D = 2 * nsites
    Q = zeros(ComplexF64, D, D)

    for ny in 0:(Ly - 1), nx in 0:(Lx - 1)
        kx = 2π * nx / Lx
        ky = 2π * ny / Ly

        Qk = qwz_bloch_grading(kx, ky; A=A, B=B, m=m)

        for y1 in 1:Ly, x1 in 1:Lx
            i0 = 2 * ((y1 - 1) * Lx + (x1 - 1))

            for y2 in 1:Ly, x2 in 1:Lx
                j0 = 2 * ((y2 - 1) * Lx + (x2 - 1))

                phase = exp(
                    im * (
                        kx * (x1 - x2) +
                        ky * (y1 - y2)
                    )
                ) / nsites

                @views Q[
                    (i0 + 1):(i0 + 2),
                    (j0 + 1):(j0 + 2)
                ] .+= phase .* Qk
            end
        end
    end

    # Remove accumulated floating-point non-Hermiticity.
    Q = (Q + Q') / 2

    return Q
end
```

This direct implementation is computationally expensive but conceptually clear. It should be used first on modest lattices. An FFT-based or displacement-block construction can be introduced after the reference version has been validated.

Essential checks are

```julia
hermiticity_error = opnorm(Q - Q')
involution_error  = opnorm(Q * Q - I)
```

For a complete periodic representation, both should be close to machine precision. For a grading constructed in a larger periodic system and then spatially truncated, (Q^2=I) will acquire boundary errors; these should decrease as the localiser centre is moved farther from the boundary and the box size grows.

---

## 11. Alternative real-space construction in the clean model

Because the indirect perturbation is known explicitly, one may construct the auxiliary direct-gap Hamiltonian (h_0) by removing its identity-channel contribution:

[
h_0=H_\gamma-g_\gamma.
]

Then diagonalise

[
h_0=U,\operatorname{diag}(\varepsilon_j),U^\dagger
]

and define

[
Q_0
===

U,\operatorname{diag}(\operatorname{sgn}\varepsilon_j),U^\dagger.
]

```julia
function grading_from_gapped_auxiliary(
    h::AbstractMatrix;
    gap_tol::Real = 1e-10
)
    F = eigen(Hermitian(Matrix(h)))

    minimum(abs, F.values) > gap_tol ||
        error("Auxiliary Hamiltonian is not gapped at zero.")

    qvals = map(F.values) do λ
        λ > 0 ? 1.0 : -1.0
    end

    Q = F.vectors * Diagonal(qvals) * F.vectors'
    Q = (Q + Q') / 2

    return Q
end
```

This construction should be compared with the Bloch/Fourier construction. They will not agree exactly near open boundaries, but they should agree increasingly well in the bulk as the system grows.

---

## 12. Compute the ordinary and graded indices side by side

For every parameter point, build the same physical Hamiltonian (H_\gamma), the same position operators, and the same localiser centre.

```julia
X, Y = build_position_operators(Lx, Ly)

Hγ = real_space_perturbed_disordered_hamiltonian_qwz(
    Lx,
    Ly;
    A=A,
    B=B,
    m=m,
    gamma=gamma,
    perturbation_type=:symmetric,
    disorder_type=:none,
    W=0.0,
    periodic_x=false,
    periodic_y=false,
    sparse_output=false
)

Q0 = clean_grading_real_space(
    Lx,
    Ly;
    A=A,
    B=B,
    m=m
)

ordinary = localiser_signature_and_gap(
    Hγ - E * I,
    X,
    Y,
    x0,
    y0;
    kappa=kappa
)

graded = localiser_signature_and_gap(
    Q0,
    X,
    Y,
    x0,
    y0;
    kappa=kappa
)
```

The resulting data row should contain at least

```julia
(
    gamma = gamma,
    E = E,
    kappa = kappa,

    ordinary_signature = ordinary.signature,
    ordinary_index = ordinary.index,
    ordinary_gap = ordinary.gap,

    graded_signature = graded.signature,
    graded_index = graded.index,
    graded_gap = graded.gap
)
```

The clean headline experiment is a (\gamma)-sweep at fixed topological (m):

* the ordinary localiser gap should decrease or close after indirect overlap;
* the graded localiser index should remain constant;
* the graded localiser gap should remain essentially independent of (\gamma);
* (Q_0) itself should be numerically identical for every (\gamma).

---

## 13. Measure the admissible (\kappa) scale

The infinite-volume invertibility estimate uses

[
M_Q=\max{|[X,Q]|,|[Y,Q]|}.
]

Numerically:

```julia
function grading_commutator_norms(Q, X, Y)
    commX = X * Q - Q * X
    commY = Y * Q - Q * Y

    MX = opnorm(commX)
    MY = opnorm(commY)

    return (
        MX = MX,
        MY = MY,
        M = max(MX, MY)
    )
end
```

The elementary sufficient upper bound is

[
\kappa<\frac{1}{2M_Q}.
]

For conservative finite-volume scans, the implementation notes propose testing an empirical window of the form

[
\frac{4}{\rho}
\lesssim
\kappa
\lesssim
\frac{1}{8M_\rho},
]

where (\rho) is the distance scale from the localiser centre to the boundary and (M_\rho) is measured on the finite system.

The lower restriction reflects finite-volume localisation: if (\kappa) is too small, the localiser does not sufficiently isolate the chosen spatial region before encountering the boundary.

The upper restriction reflects noncommutativity: if (\kappa) is too large, the position-grading commutator terms can overcome the unit grading gap.

One should not select a single (\kappa) without checking a plateau. Instead:

```julia
for kappa in kappas
    result = localiser_signature_and_gap(
        Q, X, Y, x0, y0; kappa=kappa
    )

    # Store signature and gap.
end
```

A reliable numerical value is indicated by a region in which:

1. the half-signature is constant;
2. the localiser gap remains comfortably above the numerical tolerance;
3. the value persists as the system size is increased.

---

## 14. Weak-disorder adapted-grading loop

For generic weak disorder, the target is a Hermitian involution (Q) close to (Q_0) that approximately or exactly commutes with (H).

For any current grading (Q_n), decompose the Hamiltonian into its diagonal and off-diagonal parts relative to that grading:

[
\mathcal D_{Q_n}(H)
===================

\frac12(H+Q_nHQ_n),
]

[
\mathcal O_{Q_n}(H)
===================

\frac12(H-Q_nHQ_n).
]

The residual

[
r_n
===

# |\mathcal O_{Q_n}(H)|

\frac12|[Q_n,H]|
]

measures the failure of the current grading to adapt to the disordered Hamiltonian.

A practical high-level iteration is:

```text
Input:
    disordered Hamiltonian H
    clean grading Q0
    auxiliary clean direct-gap Hamiltonian h0
    tolerance tol
    maximum iteration count maxiter

Set:
    Q = Q0

Repeat:
    1. Compute the off-diagonal residual
           R = 0.5 * (H - Q*H*Q)

    2. Express R in the positive/negative subspaces of Q0
       or of the current Q.

    3. Solve the Sylvester equation for an anti-Hermitian
       generator S:
           S_{+-} h_{--} - h_{++} S_{+-} = -R_{+-}

    4. Form the unitary update
           U = exp(S)

    5. Rotate the grading
           Q_new = U * Q * U'

    6. Restore numerical Hermiticity
           Q_new = (Q_new + Q_new') / 2

    7. Optionally purify:
           diagonalise Q_new and replace its eigenvalues
           by ±1

    8. Measure:
           residual      = 0.5 * opnorm(Q_new*H - H*Q_new)
           homotopy_dist = opnorm(Q_new - Q0)
           invol_error   = opnorm(Q_new^2 - I)

    9. Accept convergence when residual < tol.

    10. Abort or flag loss of certification if
            homotopy_dist >= 1,
        the Sylvester solve becomes ill-conditioned,
        or the residual ceases to contract.

Output:
    adapted grading Q
    convergence history
    homotopy certificate
```

In Julia-like pseudocode:

```julia
function adapted_grading_fixed_point(
    H,
    Q0,
    h0;
    tol=1e-8,
    maxiter=100
)
    Q = copy(Q0)
    history = NamedTuple[]

    for iteration in 1:maxiter
        residual_op = 0.5 * (H - Q * H * Q)
        residual = opnorm(residual_op)

        homotopy_dist = opnorm(Q - Q0)
        involution_error = opnorm(Q * Q - I)

        push!(history, (
            iteration=iteration,
            residual=residual,
            homotopy_dist=homotopy_dist,
            involution_error=involution_error
        ))

        residual < tol &&
            return (
                Q=Q,
                converged=true,
                certified=homotopy_dist < 1,
                history=history
            )

        homotopy_dist < 1 ||
            return (
                Q=Q,
                converged=false,
                certified=false,
                history=history
            )

        # Project residual into opposite clean-band sectors.
        Pminus = (I - Q0) / 2
        Pplus  = (I + Q0) / 2

        hminus = Pminus * h0 * Pminus
        hplus  = Pplus  * h0 * Pplus
        target = Pplus * residual_op * Pminus

        # Conceptual step:
        # Solve S_pm*hminus - hplus*S_pm = -target.
        S_pm = solve_band_sylvester(
            hplus,
            hminus,
            -target
        )

        S = S_pm - S_pm'
        U = exp(S)

        Q = U * Q * U'
        Q = (Q + Q') / 2

        # Optional spectral purification.
        F = eigen(Hermitian(Q))
        signs = map(λ -> λ >= 0 ? 1.0 : -1.0, F.values)
        Q = F.vectors * Diagonal(signs) * F.vectors'
    end

    return (
        Q=Q,
        converged=false,
        certified=false,
        history=history
    )
end
```

This is intentionally schematic: the Sylvester solve must be carried out in explicit positive- and negative-band bases, rather than on matrices padded with null sectors. The implementation plan proposes developing both the proof-faithful (S)-iteration and a more pragmatic self-consistent projector iteration, then checking their agreement at small disorder.

---

## 15. A simpler self-consistent projector prototype

Before implementing the full proof-faithful rotation, a useful exploratory algorithm is:

1. start with (Q_0);
2. take the part of (H) diagonal in the current grading;
3. remove any identified scalar or trace component;
4. flatten this auxiliary direct-gap operator;
5. repeat.

```julia
function self_consistent_grading(
    H,
    Q0;
    tol=1e-8,
    maxiter=100,
    gap_tol=1e-8
)
    Q = copy(Q0)
    history = NamedTuple[]

    for iteration in 1:maxiter
        # Part of H diagonal with respect to Q.
        h = 0.5 * (H + Q * H * Q)

        # In a two-orbital translationally structured calculation,
        # remove the identified identity-channel component here.
        h_traceless = remove_scalar_channel(h)

        F = eigen(Hermitian(h_traceless))

        minimum(abs, F.values) > gap_tol ||
            return (
                Q=Q,
                converged=false,
                reason=:direct_gap_closed,
                history=history
            )

        signs = map(λ -> λ > 0 ? 1.0 : -1.0, F.values)
        Qnew = F.vectors * Diagonal(signs) * F.vectors'
        Qnew = (Qnew + Qnew') / 2

        update = opnorm(Qnew - Q)
        residual = 0.5 * opnorm(Qnew * H - H * Qnew)
        homotopy_dist = opnorm(Qnew - Q0)

        push!(history, (
            iteration=iteration,
            update=update,
            residual=residual,
            homotopy_dist=homotopy_dist
        ))

        Q = Qnew

        update < tol &&
            return (
                Q=Q,
                converged=true,
                certified=homotopy_dist < 1,
                history=history
            )
    end

    return (
        Q=Q,
        converged=false,
        certified=false,
        history=history
    )
end
```

This prototype is not a substitute for the mathematically specified fixed-point map unless their equivalence is demonstrated. Its role is to:

* expose programming and basis-management errors;
* provide an independent small-(W) comparison;
* test whether a stable adapted grading appears numerically;
* identify the disorder scale at which self-consistency breaks down.

---

## 16. Diagnostics that must accompany every graded result

A half-signature should never be stored by itself. For each computed grading (Q), record:

[
\epsilon_{\mathrm{herm}}
========================

|Q-Q^\dagger|,
]

[
\epsilon_{\mathrm{inv}}
=======================

|Q^2-I|,
]

[
\epsilon_{\mathrm{comm}}
========================

\frac12|[Q,H]|,
]

[
d_{\mathrm{hom}}
================

|Q-Q_0|,
]

[
M_\rho
======

\max{|[X,Q]|,|[Y,Q]|},
]

[
g_{\mathrm{loc}}
================

\min_{\lambda\in\sigma(\widetilde L_\kappa)}|\lambda|.
]

A result in the weak-disorder calculation should be labelled, for example,

```julia
certified = (
    converged &&
    homotopy_distance < 1 &&
    involution_error < involution_tol &&
    localiser_gap > localiser_gap_tol &&
    kappa >= kappa_lower &&
    kappa <= kappa_upper
)
```

The word “certified” here means certified as belonging to the intended numerical construction and homotopy class. It should not be used to imply that the still-pending disordered finite-volume signature theorem has already been proved.

---

## 17. Disorder realisations must be held fixed during parameter comparisons

The current parameter-sweep code generates a fresh disorder realisation each time the Hamiltonian builder is called. That prevents meaningful comparisons between different (\kappa), (\gamma), position, or disorder strengths.

Instead, generate a dimensionless disorder profile once:

```julia
using Random

function disorder_profile(
    Lx::Int,
    Ly::Int;
    seed::Int
)
    rng = MersenneTwister(seed)
    return rand(rng, Lx, Ly) .- 0.5
end
```

Then define

[
V(W)=W,V_{\mathrm{shape}}.
]

```julia
profile = disorder_profile(Lx, Ly; seed=seed)

for W in Ws
    onsite_disorder = W .* profile
    H = build_hamiltonian_with_given_disorder(
        ...,
        onsite_disorder=onsite_disorder
    )

    # All κ and position scans at this W use this same H.
end
```

This produces a continuous disorder path in (W) and allows:

* a meaningful convergence threshold to be located;
* a meaningful homotopy-distance curve to be plotted;
* localiser gaps at different (\kappa) to be compared for the same system;
* seed-to-seed variation to be separated from parameter dependence.

This refactoring is a correctness requirement, not merely an optimisation.

---

## 18. Minimal validation sequence

The implementation should be accepted in the following order.

### Test A: algebraic validity of the clean grading

Verify

[
Q^\dagger\approx Q,
\qquad
Q^2\approx I.
]

Verify that (Q) is independent of (\gamma).

### Test B: clean topological benchmark

Compute a momentum-space Chern number independently, for example using the Fukui–Hatsugai–Suzuki lattice formula.

Check

[
\frac12\operatorname{sig}\widetilde L_\kappa
============================================

C_{\mathrm{FHS}}.
]

### Test C: ordinary-versus-graded (\gamma)-sweep

At fixed topological (m), increase (\gamma) from the insulating regime into deep indirect overlap.

Expected behaviour:

[
\operatorname{gap}L_\kappa(H_\gamma-E)
\longrightarrow 0
]

for the ordinary localiser, while

[
\operatorname{gap}\widetilde L_\kappa(Q_0)
]

and

[
\frac12\operatorname{sig}\widetilde L_\kappa(Q_0)
]

remain constant up to finite-size and numerical errors.

### Test D: (\kappa)-plateau

Scan (\kappa) logarithmically. Identify an interval on which:

[
\frac12\operatorname{sig}\widetilde L_\kappa=C
]

and the localiser gap is nonzero.

Compare the upper edge with the measured commutator scale (1/M_\rho).

### Test E: spatial and size stability

Move ((x_0,y_0)) from the centre toward the boundary and repeat for increasing (L_x,L_y).

The interior index should stabilise first. Boundary effects should retreat outward as the system size grows.

### Test F: weak-disorder grading convergence

For fixed disorder shape, increase (W) and record:

[
|\mathcal O_Q(H)|,
\qquad
|Q-Q_0|,
\qquad
|Q^2-I|,
\qquad
\operatorname{gap}\widetilde L_\kappa.
]

Within the weak-disorder convergence regime, the expected numerical pattern is:

* residual decreases under iteration;
* (|Q-Q_0|<1);
* the graded half-signature remains on the clean topological value;
* the graded localiser remains gapped over a measurable (\kappa)-plateau.

At larger (W), the iteration may stop contracting, the homotopy distance may approach one, or the auxiliary direct gap may become ill-conditioned. That point should be reported as the empirical boundary of the present construction, not automatically as a topological transition.

---

## 19. Interpretation of the output

The final number is computed in exactly the same way as for the ordinary spectral localiser:

[
\boxed{
C_{\mathrm{numerical}}
======================

\frac12
\left(
N_+(\widetilde L_\kappa)
------------------------

N_-(\widetilde L_\kappa)
\right).
}
]

What has changed is the object supplied to the localiser.

For the ordinary method,

[
C_{\mathrm{ordinary}}
=====================

\frac12\operatorname{sig}L_\kappa(H-E),
]

so the result is meaningful when energy separates the states of interest.

For the new method,

[
C_{\mathrm{graded}}
===================

\frac12\operatorname{sig}\widetilde L_\kappa(Q),
]

so the result is meaningful when a local band grading separates the states of interest.

The conceptual progression is therefore

[
\text{physical Hamiltonian}
\quad\longrightarrow\quad
\text{band projector}
\quad\longrightarrow\quad
\text{local grading}
\quad\longrightarrow\quad
\text{graded localiser}
\quad\longrightarrow\quad
\text{half-signature}.
]

The ordinary localiser asks whether topology is visible through a gap in the **energy spectrum**. The repaired localiser asks whether topology is visible through a gap in a local **band-label operator**. In an indirect-gap system, the latter is the appropriate question.
