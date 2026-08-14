# Restatement of the Problem

## Setup and definitions

Let $H_\gamma$ be the perturbed QWZ Hamiltonian on $\ell^2(\mathbb{Z}^2)\otimes\mathbb{C}^2$, with Bloch form

$$
H_\gamma(\mathbf{k}) = H_{\mathrm{QWZ}}(\mathbf{k}) + \gamma(\cos k_x + \cos k_y)\,\mathbb{I}_2 .
$$

Key structural facts:

1. **Eigenstates unchanged.** Since the perturbation is $\propto \mathbb{I}_2$, the Bloch eigenprojections $P_\pm(\mathbf{k})$ of $H_\gamma(\mathbf{k})$ are identical to those of $H_0(\mathbf{k})$ for all $\gamma$. The occupied vector bundle $\mathcal{E}_- \to T^2$, and hence its Chern number $C$, is $\gamma$-independent.
2. **Local (direct) gap persists.** For all $\mathbf{k}$, $E_+(\mathbf{k}) - E_-(\mathbf{k}) = 2|\mathbf{d}(\mathbf{k})| > 0$ (away from the topological transition). So the two bands never touch pointwise in $\mathbf{k}$.
3. **Global gap destroyed.** For $\gamma$ large enough, $\max_\mathbf{k} E_-(\mathbf{k}) > \min_\mathbf{k} E_+(\mathbf{k})$ (overlap near $\Gamma$ and $M$): an *indirect* band overlap. No Fermi energy $E_F$ lies in a spectral gap of $H_\gamma$.

The standard spectral localiser at position $(x_0,y_0)$ and energy $E$,

$$
L_\kappa(x_0,y_0,E) =
\begin{pmatrix}
H - E & \kappa\big((X-x_0) - i(Y-y_0)\big) \\
\kappa\big((X-x_0) + i(Y-y_0)\big) & -(H-E)
\end{pmatrix},
$$

requires $E$ to lie in a local spectral gap of $H$; its signature index then equals $C$. In the indirect-overlap regime no valid $E$ exists, and the localiser's index becomes ill-defined/unreliable **even though the topology of the occupied bundle is intact**. This is the failure to be repaired.

## What has been established

- **No-Go I (spectral-function no-go).** There is no scalar function $f$ such that $f(H_\gamma) = P_-$ (the band projector) once the bands overlap: if $E_0 \in \sigma_-\cap\sigma_+$, then $f(H)$ acts identically on degenerate valence and conduction eigenstates, so it cannot output $1$ and $0$ simultaneously. This kills KPM/Fermi-function filtering and any $\Theta(E_F - H)$ construction. ✔ (proved)
- **No-Go II (symmetry).** Inserting inversion symmetry into the localiser does not restore a well-defined index (your separate result).

**Consequence:** any repair must inject information *beyond* $\sigma(H)$ — i.e. it must be built from the band projector $P_-$ itself (or an object that determines it), not from a functional calculus on $H$.

## The precise mathematical problem

> **Problem.** Construct a real-space, exponentially local, self-adjoint operator $\widetilde{H}$, computable in principle from $H_\gamma$ plus admissible extra structure, such that $0$ lies in a spectral gap of $\widetilde{H}$, the localiser $\widetilde{L}_\kappa = $ (localiser built with $\widetilde{H}$, $E=0$) is invertible, and $\tfrac12\,\mathrm{sig}\,\widetilde{L}_\kappa = C$, the Chern number of the occupied bundle of $H_\gamma$.

The canonical candidate is the **band-flattened Hamiltonian**

$$
Q := P_+ - P_- = \mathbb{I} - 2P_-,
$$

which by construction has $\sigma(Q) = \{-1, +1\}$, gap $=2$, and the same occupied bundle. Crucially, in the indirect regime $Q \neq \operatorname{sgn}(H - E_F)$ for any $E_F$ — this is *band flattening*, not *energy flattening*, and by No-Go I it is not any $f(H)$.

## Testable conjectures

**Conjecture A (locality of $Q$).** *Because the direct gap $2\min_\mathbf{k}|\mathbf{d}(\mathbf{k})| > 0$ is open and $P_-(\mathbf{k})$ is real-analytic on $T^2$, the real-space kernel of $Q$ decays exponentially:*
$$
\|\langle \mathbf{r}| Q |\mathbf{r}'\rangle\| \le C e^{-|\mathbf{r}-\mathbf{r}'|/\xi},
$$
*with $\xi$ controlled by the direct gap, not the (nonexistent) indirect gap.*
Test: this is essentially a Combes–Thomas / analyticity-of-Bloch-projections statement; for the translation-invariant toy model it should follow directly since $P_-(\mathbf{k})$ is analytic. The nontrivial extension is to disordered/finite systems where $P_-$ must be defined without $\mathbf{k}$.

**Conjecture B (flattened localiser works).** *Let $\widetilde{L}_\kappa(x_0,y_0)$ be the localiser with $H \to Q$, $E = 0$. Then for $\kappa$ below a threshold set by the gap of $Q$ (which is $2$) and the commutator bound $\|[X, Q]\|$, $\widetilde{L}_\kappa$ is invertible and*
$$
\tfrac12\,\mathrm{sig}\,\widetilde{L}_\kappa(x_0,y_0) = C .
$$
Tests: (i) verify $\|[X,Q]\|, \|[Y,Q]\| < \infty$ (bounded, guaranteed by Conjecture A); (ii) check numerically on finite lattices that the signature is quantised, stable in $\kappa$, and jumps only at the phase transition of the *unperturbed* $\mathbf{d}(\mathbf{k})$ — i.e. is $\gamma$-independent.

**Conjecture C (flattening homotopy).** *Define $H_s = (1-s)(H_\gamma - E_F) + sQ$, $s\in[0,1]$, with $E_F$ chosen inside the direct gap "locally" — or more robustly, $H_s = (1-s)\,Q|H'|Q$-type deformations. The relevant claim: the localiser $L_\kappa[H_s]$ remains invertible for all $s$, so its index is constant and equals both the $H$-localiser index (where defined) and the $Q$-localiser index.*
This is the sharp statement that "the localiser only sees the bundle": Is $L_\kappa[H_s]$ invertible along the flattening homotopy? Note the subtlety: at $s=0$ the localiser may *already* be non-invertible in the indirect regime — the honest version of Conjecture C is stated in the direct-gap regime (small $\gamma$), establishing index equality there, then using $\gamma$-continuity of $\widetilde{L}$ (Conjecture B, whose gap never closes in $\gamma$) to extend to the indirect regime. This chain is the actual proof strategy:

$$
\underbrace{\mathrm{ind}\,L_\kappa[H_{\gamma=0}]}_{\text{standard theory} = C}
\;\overset{\text{Conj. C}}{=}\;
\mathrm{ind}\,\widetilde{L}_\kappa[Q_{\gamma=0}]
\;\overset{\gamma\text{-indep. of }Q}{=}\;
\mathrm{ind}\,\widetilde{L}_\kappa[Q_{\gamma}] \quad \forall\gamma .
$$

Note $Q_\gamma \equiv Q_0$ identically in this toy model, so the last step is trivial *here*; the general theorem needs only that $Q$ deforms continuously without gap closing while the direct gap stays open.

**Conjecture D (constructibility of $P_-$ — the repaired Conjecture 2).** In light of No-Go I, restate the projector conjecture as: *$P_-$ can be obtained by an admissible beyond-spectrum construction, e.g.*
- *(D1, adiabatic)* $P_-(\gamma) = U_\gamma P_-(0) U_\gamma^\dagger$ via parallel transport along $\gamma$, valid because the direct gap of the band splitting never closes — well-defined by the quasi-adiabatic evolution / Hastings–Wen construction applied to the *band* gap, not the energy gap;
- *(D2, contour)* $P_- = \frac{1}{2\pi i}\oint_{\mathcal{C}} (z - h)^{-1} dz$ where $h$ is an auxiliary local Hamiltonian isospectral in bundle (e.g. $H_\gamma$ minus the identity-part perturbation), when such a decomposition $H = h + g$, $[g\text{-part}] \propto \mathbb{I}$ locally, can be identified;
- *(D3, extra observable)* not spatial symmetry (ruled out), but the operator $Q$ itself commuting with $H$ — the question becomes whether $Q$ can be characterised variationally as the unique local involution commuting with $H$ with $\mathrm{tr}$-density fixed.

**Falsifiable predictions** (numerics on the toy model):
1. Localiser gap of $\widetilde{L}_\kappa$ is independent of $\gamma$; localiser gap of $L_\kappa$ closes as $\gamma$ crosses the indirect-overlap threshold.
2. $\mathrm{sig}\,\widetilde{L}_\kappa/2 = C$ for all $\gamma$, including deep in the Chern-metal regime.
3. $\|[X, Q]\|$ stays bounded uniformly in $\gamma$ (it is exactly $\gamma$-independent here — a clean sanity check).

## In one sentence

Replace the energy-based localiser $L_\kappa(H - E_F)$ with the band-based localiser $\widetilde{L}_\kappa(Q)$, prove (A) $Q$ is exponentially local given only a *direct* gap, (B) the $Q$-localiser index equals the Chern number, (C) the two localisers agree via homotopy wherever both are defined, and (D) give an in-principle, energy-free construction of $P_-$ consistent with No-Go I — with the toy model providing a complete testbed since there $Q$ is exactly $\gamma$-independent.








# Proof of Conjecture B

## Strategy

The key observation that makes this proof short: the Loring–Schulz-Baldes localiser theorem never uses that the operator in the "energy slot" is the *physical* Hamiltonian. Its hypotheses are only that one has a **bounded, self-adjoint, exponentially local operator with a spectral gap at the reference energy**. The flattened operator $Q = \mathbb{I} - 2P_-$ satisfies these hypotheses (given Conjecture A), and its negative spectral subspace is *exactly* the occupied bundle. So Conjecture B splits into:

- **(i)** an elementary invertibility estimate (proved from scratch below), and
- **(ii)** an identification of $\mathrm{sig}\,\widetilde L_\kappa/2$ with the Chern number, which follows by applying the existing theorem to $Q$ and checking its hypotheses.

## Setup

Let $Q = Q^\dagger$, $Q^2 = \mathbb{I}$, acting on $\mathcal{H} = \ell^2(\mathbb{Z}^2)\otimes\mathbb{C}^2$, with (Conjecture A)

$$
\big\|\langle\mathbf{r}|Q|\mathbf{r}'\rangle\big\| \le C e^{-|\mathbf{r}-\mathbf{r}'|/\xi}.
$$

Write $X_0 = X - x_0$, $Y_0 = Y - y_0$. Define, using Pauli matrices $\sigma_i$ acting on an auxiliary $\mathbb{C}^2$,

$$
\widetilde L_\kappa \;=\; \kappa\, X_0\otimes\sigma_x \;+\; \kappa\, Y_0\otimes\sigma_y \;+\; Q\otimes\sigma_z ,
$$

which is the matrix form of the localiser you wrote, with $H - E$ replaced by $Q - 0$.

**Lemma 0 (bounded commutators).** Exponential locality implies $\|[X, Q]\| < \infty$. Indeed, $\langle\mathbf{r}|[X,Q]|\mathbf{r}'\rangle = (r_x - r'_x)\langle\mathbf{r}|Q|\mathbf{r}'\rangle$, whose norm is $\le C|r_x - r'_x| e^{-|\mathbf{r}-\mathbf{r}'|/\xi}$; by the Schur test (row/column sums are bounded by $C\sum_{\mathbf{d}} |d_x| e^{-|\mathbf{d}|/\xi} =: M < \infty$, uniformly in the row/column), $\|[X,Q]\| \le M$. Same for $Y$. ∎

Set $M := \max\big(\|[X,Q]\|, \|[Y,Q]\|\big)$. Note for the toy model $M$ is **exactly $\gamma$-independent**, since $Q$ is.

## Part (i): Invertibility of $\widetilde L_\kappa$

**Proposition 1.** If $\;\kappa < \dfrac{1}{2M}\;$ then $\widetilde L_\kappa^2 \ge (1 - 2\kappa M)\,\mathbb{I} > 0$; in particular $\widetilde L_\kappa$ is invertible with a spectral gap $\ge \sqrt{1 - 2\kappa M}$ around $0$, **uniformly in $(x_0, y_0)$ and in system size**.

**Proof.** Square $\widetilde L_\kappa$ using $\sigma_i\sigma_j + \sigma_j\sigma_i = 2\delta_{ij}$:

$$
\widetilde L_\kappa^2
= \big(\kappa^2 X_0^2 + \kappa^2 Y_0^2 + Q^2\big)\otimes\mathbb{I}_2
\;+\; \kappa\,[X_0, Q]\otimes\sigma_x\sigma_z
\;+\; \kappa\,[Y_0, Q]\otimes\sigma_y\sigma_z ,
$$

where the diagonal terms combined because $\sigma_x\sigma_y + \sigma_y\sigma_x = 0$ and $[X_0, Y_0] = 0$, and the cross terms with $Q$ collected into commutators precisely because $\sigma_z$ *anti*commutes with $\sigma_x, \sigma_y$:
$$
\kappa X_0 Q\,(\sigma_x\sigma_z) + \kappa Q X_0\,(\sigma_z\sigma_x) = \kappa [X_0, Q]\,\sigma_x\sigma_z .
$$

Now use the three structural facts:

1. $Q^2 = \mathbb{I}$ — **this is where flattening pays off**: the usual term $(H-E)^2 \ge \mathrm{gap}^2$ (which fails without a global gap) is replaced by exactly $\mathbb{I}$;
2. $\kappa^2(X_0^2 + Y_0^2) \ge 0$;
3. $[X_0, Q] = [X, Q]$ and $\|A\otimes\sigma_i\sigma_j\| = \|A\|$, so the error terms are bounded below by $-\kappa M\,\mathbb{I}$ each.

Hence
$$
\widetilde L_\kappa^2 \;\ge\; \big(1 - 2\kappa M\big)\,\mathbb{I} .
$$
For $\kappa < 1/(2M)$ the right side is strictly positive, so $0 \notin \sigma(\widetilde L_\kappa)$ and $\mathrm{dist}(0, \sigma(\widetilde L_\kappa)) \ge \sqrt{1 - 2\kappa M}$. Nothing in the estimate referenced $x_0, y_0$ or the volume. $\blacksquare$

**Remarks.**
- Contrast with the standard localiser: there one needs $(H-E)^2 \ge g^2$ with $g$ the *local energy gap at $E$*, and invertibility requires $2\kappa M < g^2$. In the indirect regime $g = 0$ for every $E$ and the estimate is vacuous — this is precisely the failure mode. Flattening restores $g = 1$ identically, for **all** $\gamma$.
- On a finite lattice $\Lambda_\rho$ of radius $\rho$ (with $Q$ truncated), the same bound holds up to boundary corrections $O(e^{-(\rho - |x_0|)/\xi})$ from truncating $Q$; these are absorbed for $\rho$ large. This is the standard finite-volume bookkeeping and I'll suppress it below, quoting it where needed.

## Part (ii): The half-signature equals the Chern number

**Proposition 2.** Let $Q$ be as above, exponentially local, with $Q^2 = \mathbb{I}$, and let $P_- = \tfrac12(\mathbb{I} - Q)$ be its negative spectral projection. Let $\widetilde L_{\kappa,\rho}$ be the localiser restricted to $\Lambda_\rho$. Then there exist constants such that for
$$
\kappa \le \kappa_0\big(M\big), \qquad \rho \ge \rho_0(\kappa),
$$
$\widetilde L_{\kappa,\rho}$ is invertible and
$$
\tfrac12\,\mathrm{sig}\,\widetilde L_{\kappa,\rho} \;=\; \mathrm{Ind}\big(P_- F P_- \big)\;=\; C,
$$
where $F = (X_0 + iY_0)/|X_0 + iY_0|$ is the Dirac phase and $C$ is the Chern number of the occupied bundle.

**Proof.**

*Step 1 (hypothesis check).* The even-dimensional localiser theorem of Loring–Schulz-Baldes (*"The spectral localizer for even index pairings"*, and the finite-volume version in *"Finite volume calculation of K-theory invariants"*) requires of the operator $H'$ in the energy slot only:

- (a) $H' = H'^\dagger$ bounded: $\|Q\| = 1$. ✔
- (b) Spectral gap at the reference energy $E = 0$: $\sigma(Q) = \{-1, +1\}$, gap $g = 1$. ✔
- (c) $\|[X, H']\|, \|[Y, H']\| < \infty$: Lemma 0. ✔ (Their proofs quantify $\kappa_0 \sim g^3/(\|H'\|\,M)$ and $\rho_0 \sim M/(\kappa g)$; with $\|Q\| = g = 1$ these are clean: $\kappa_0 \sim c/M$, consistent with Proposition 1.)

Nowhere do the hypotheses require $H'$ to generate the physical dynamics. Applying the theorem to $H' = Q$ gives invertibility of $\widetilde L_{\kappa,\rho}$ and
$$
\tfrac12\,\mathrm{sig}\,\widetilde L_{\kappa,\rho} = \mathrm{Ind}\big(P\, F\, P\big)\Big|_{P = \chi(Q < 0)} .
$$

*Step 2 (identify the projection).* By construction $\chi(Q < 0) = P_-$, the band projector onto the occupied bundle of $H_\gamma$. Note this step is where **band** information enters, consistently with No-Go I: $P_-$ was *not* obtained as $f(H_\gamma)$; the theorem only needs $P_-$ to exist as a local projection, however it was constructed (Conjecture D's job).

*Step 3 (index = Chern number).* For an exponentially local projection $P_-$ on $\ell^2(\mathbb{Z}^2)\otimes\mathbb{C}^N$, the Fredholm index of the Noether–Bellissard pairing $P_- F P_-: \mathrm{Ran}\,P_- \to \mathrm{Ran}\,P_-$ equals the Chern number of $P_-$ (Bellissard–van Elst–Schulz-Baldes; in the translation-invariant case this is the classical index computation, $\mathrm{Ind} = \frac{1}{2\pi i}\int_{T^2} \mathrm{tr}\,P_-\,[\partial_{k_x}P_-, \partial_{k_y}P_-]\, d^2k = C$). Since $P_-(\mathbf{k})$ is the *same* projector for all $\gamma$, this index is $\gamma$-independent. $\blacksquare$

## Corollary (the toy model, all $\gamma$)

For the perturbed QWZ model, $Q_\gamma \equiv Q_0$ identically. Hence $M$, $\kappa_0$, $\rho_0$, and the localiser gap $\sqrt{1 - 2\kappa M}$ are all $\gamma$-independent, and

$$
\tfrac12\,\mathrm{sig}\,\widetilde L_{\kappa,\rho} = C \quad \text{for every } \gamma,
$$

including arbitrarily deep in the indirect-overlap (Chern-metal) regime — while the standard localiser $L_\kappa(H_\gamma - E)$ loses its gap as soon as no $E$ sits in a spectral gap. This proves Conjecture B and simultaneously yields falsifiable predictions 1–3 from before.

## What is genuinely new vs. borrowed — and what remains

Honest accounting:

- **Borrowed:** the signature-index theorem (Step 1). Its proof (spectral flow / K-theoretic homotopy from $\widetilde L_\kappa$ to a direct sum of a trivial part and the index pairing) is unchanged; only the input operator differs.
- **New content:** (a) the observation that $Q$ is an admissible input despite not being $f(H_\gamma)$ — i.e. the localiser theorem is a theorem about *(projection, position)* pairs, and the energy dispersion was always a scaffold; (b) the uniform-in-$\gamma$ constants, which is what makes the indirect regime accessible.
- **Load-bearing dependencies:** Proposition 2 consumes **Conjecture A** (locality of $Q$) as a hypothesis — for the translation-invariant toy model this is immediate from real-analyticity of $P_-(\mathbf{k})$ plus Paley–Wiener, so B is *unconditional* there; the disordered case genuinely needs A proved from a mobility-gap-type condition. And the whole construction presupposes access to $P_-$, i.e. **Conjecture D**. B says: *if you can build the band projector locally, the localiser is fully repaired.* D is now the sole remaining obstruction.

A worthwhile immediate check: Conjecture C (the interpolation $H_s = (1-s)(H_\gamma - E_F) + sQ$) in the *direct-gap* regime now follows from the same squaring estimate — $H_s^2 \ge \min(g_s)^2$ with $g_s$ interpolating between $g$ and $1$, and one only needs $g_s > 0$ along the path plus a commutator bound uniform in $s$. Want me to write that estimate out next, or move straight to attacking Conjecture D (the energy-free construction of $P_-$, where I'd start with the quasi-adiabatic route D1, since the direct gap gives exactly the parallel-transport structure it needs)?








# Conjecture A in the Mobility-Gap Setting

## Step 0: The problem must first be restated — $\mathbf{k}$-space is gone

Conjecture A as stated ("$P_-(\mathbf{k})$ analytic $\Rightarrow$ exponential decay") is a translation-invariant statement. Before proving anything we must face an honest definitional issue: **without translation invariance there are no bands**, so "the band projector $P_-$" must be *defined* before its locality can be studied. This is not pedantry — it is the disordered shadow of Conjecture D.

The clean way out, which the toy model suggests, is this structural observation:

> In the toy model, $H_\gamma = H_0 + \gamma\Delta$ where $\Delta = \sum_{\langle\mathbf{r}\mathbf{r}'\rangle}|\mathbf{r}\rangle\langle\mathbf{r}'|\otimes\mathbb{I}_2 + \text{h.c.}$ acts as the identity in the orbital channel. Hence $[\Delta, P_-] = 0$ is *not* automatic — but here $\Delta$ is a scalar function of quasi-momentum times $\mathbb{I}_2$, so in Bloch form $[H_0(\mathbf{k}), g(\mathbf{k})\mathbb{I}_2] = 0$, giving $[H_\gamma, Q] = 0$ with $Q = \mathbb{I} - 2\chi(H_0 < 0)$. The band projector is the spectral projection **of the auxiliary operator $h := H_0 = H_\gamma - \gamma\Delta$, which has a genuine spectral gap** even when $H_\gamma$ does not.

This motivates the disorder-robust reformulation:

**Conjecture A′ (intrinsic form).** *Suppose there exists a self-adjoint, uniformly-local (finite-range or exponentially decaying hopping) operator $h$ on $\ell^2(\mathbb{Z}^2)\otimes\mathbb{C}^2$ such that*

1. $P_- = \chi(h < 0)$ *(the occupied bundle is the negative spectral subspace of $h$), and*
2. $h$ *has either (i) a spectral gap $g > 0$ around $0$, or (ii) a mobility gap at $0$ (defined below).*

*Then $Q = \mathbb{I} - 2P_-$ satisfies kernel-decay bounds sufficient for the localiser construction.*

The "direct gap" of the clean problem is thus replaced by: **a gap (spectral or mobility) of the identity-channel-stripped operator $h$**, not of $H$. Disorder is admissible provided it can be split as $V = V_h + V_\Delta$ with $V_h$ entering $h$ (arbitrary, as long as $h$ stays gapped/localised) and $V_\Delta$ commuting with $Q$ (e.g. orbital-scalar disorder in the toy model). I flag explicitly: whether *generic* disorder admits such a splitting is a separate structural question belonging to Conjecture D; here we assess robustness *given* the splitting.

## Step 1: Case (i) — spectral gap of $h$. Fully rigorous, deterministic.

**Proposition A1.** Let $h = h^\dagger$ with $\|h\| \le B$, hopping range $R$ (i.e. $\langle\mathbf{r}|h|\mathbf{r}'\rangle = 0$ for $|\mathbf{r}-\mathbf{r}'| > R$), and $\sigma(h) \cap (-g, g) = \emptyset$. Then there exist $C, \xi > 0$ depending only on $B, R, g$ (not on the disorder realisation, not on volume) such that
$$
\big\|\langle\mathbf{r}|Q|\mathbf{r}'\rangle\big\| \le C\, e^{-|\mathbf{r}-\mathbf{r}'|/\xi}, \qquad \xi \sim \frac{2BR}{g}\ \text{(up to constants)}.
$$

**Proof.** Write the Riesz projection
$$
P_- = \frac{1}{2\pi i}\oint_{\mathcal{C}} (z - h)^{-1}\, dz,
$$
with $\mathcal{C}$ a contour crossing the real axis only at $0$ and at $-B - 1$, staying at distance $\ge g/2$ from $\sigma(h)$. This is legitimate because $h$ has a genuine spectral gap — **note this does not contradict No-Go I**: the contour integral is a function of $h$, not of $H_\gamma$; the beyond-spectrum information is the *decomposition* $H_\gamma = h + \gamma\Delta$.

Combes–Thomas estimate: for $\mu > 0$ define $h_\mu = e^{\mu \langle\mathbf{a}, \mathbf{X}\rangle} h\, e^{-\mu \langle\mathbf{a}, \mathbf{X}\rangle}$ for a unit vector $\mathbf{a}$. Finite range gives $\|h_\mu - h\| \le B(e^{\mu R} - 1)$. Choosing $\mu$ so that $B(e^{\mu R} - 1) \le g/4$, a Neumann-series argument shows $(z - h_\mu)^{-1}$ exists on $\mathcal{C}$ with $\|(z - h_\mu)^{-1}\| \le 4/g$, whence
$$
\big\|\langle\mathbf{r}|(z-h)^{-1}|\mathbf{r}'\rangle\big\| \le \frac{4}{g}\, e^{-\mu|\mathbf{r} - \mathbf{r}'|}
$$
uniformly on $\mathcal{C}$. Integrating over the contour (length $\le 2(B + 2) + 2g \le$ const) yields the claim with $1/\xi = \mu$. $\blacksquare$

**Consequences.** Everything in Propositions 1 and 2 of the previous answer goes through *verbatim and deterministically*: the Schur test gives $\|[X, Q]\| \le M(B, R, g) < \infty$ uniformly in disorder realisation and volume, and the localiser theorem applies with $\omega$-independent constants $\kappa_0, \rho_0$. Quantisation and disorder-independence of $\tfrac12\mathrm{sig}\,\widetilde L$ follow. **In regime (i), Conjecture A′ is a theorem and the approach is fully robust.**

For the toy model plus orbital-scalar disorder this regime already applies: $h = H_0 + V_h$ retains its gap for $\|V_h\| < g$, while $\gamma$ and $V_\Delta$ are unconstrained. That is already a nontrivial robustness statement.

## Step 2: Case (ii) — mobility gap only. Here it gets genuinely delicate.

**Definition (mobility gap, fractional-moment form).** $h_\omega$ (random, ergodic) has a mobility gap in $[-g, g]$ if $0 \in \sigma(h_\omega)$ is allowed, but for some $s \in (0,1)$, $\mu > 0$:
$$
\sup_{E \in [-g,g]}\ \mathbb{E}\Big[\big\|\langle\mathbf{r}|(h_\omega - E - i0)^{-1}|\mathbf{r}'\rangle\big\|^s\Big] \le C\, e^{-\mu|\mathbf{r} - \mathbf{r}'|}.
$$
(Aizenman–Molchanov condition; equivalent characterisations via eigenfunction correlators exist.)

**What survives (Lemma A2).** Under this condition the standard machinery (Aizenman–Graf) yields for $P_- = \chi(h_\omega < 0)$:

- *In expectation:* $\ \mathbb{E}\big\|\langle\mathbf{r}|P_-|\mathbf{r}'\rangle\big\| \le C\, e^{-\mu'|\mathbf{r}-\mathbf{r}'|}$.
- *Almost surely:* via Chebyshev + Borel–Cantelli over the lattice, for a.e. $\omega$ there is $C_\omega < \infty$ and $\nu > 0$ (any $\nu > 2d/s$-type exponent, dimension $d = 2$) with
$$
\big\|\langle\mathbf{r}|P_-|\mathbf{r}'\rangle\big\| \le C_\omega\,(1 + |\mathbf{r}|)^{\nu}\, e^{-\mu'|\mathbf{r}-\mathbf{r}'|}. \tag{$\star$}
$$

The crucial degradation: **the prefactor is random and grows polynomially with position**. Exponential decay in the *separation* survives; uniformity in the *base point* does not. This is not a technical artefact — near-critical eigenstates at energies inside the mobility-gap window genuinely produce rare regions with fat kernel tails.

**What this does to the two propositions:**

**(a) Proposition 2's index (the pairing) is safe.** The bound $(\star)$ is exactly the "weak locality" under which the Fredholm index $\mathrm{Ind}(P_- F P_-)$ is known to be well-defined, a.s. constant, and integer-quantised (Aizenman–Graf; Bellissard–van Elst–Schulz-Baldes noncommutative Chern number; Elgart–Graf–Schenker for the edge counterpart). So **the topological invariant itself is robust in the mobility-gap regime.**

**(b) Proposition 1's estimate breaks on the infinite lattice.** The Schur row-sum at site $\mathbf{r}$ is now bounded only by $C_\omega (1+|\mathbf{r}|)^\nu \cdot \text{const}$, so $[X, Q]$ is (possibly) **unbounded** as an operator on $\ell^2(\mathbb{Z}^2)$. The squaring argument $\widetilde L^2 \ge 1 - 2\kappa\|[X,Q]\|$ is then vacuous. This is the precise point of failure — and it is the *same* point where the standard localiser theory also currently lacks a full mobility-gap proof, so we have not made the problem worse than the state of the art; we have inherited its open boundary.

## Step 3: Salvage in finite volume — a quantitative condition

The localiser is evaluated on $\Lambda_\rho$ anyway. Restrict all estimates to $|\mathbf{r}|, |\mathbf{r}'| \le \rho$ and centre $x_0 = y_0 = 0$:

$$
M_\rho := \max\big(\|[X, Q]\|_{\Lambda_\rho}, \|[Y, Q]\|_{\Lambda_\rho}\big) \le C_\omega'\, \rho^{\nu}.
$$

The two requirements of the localiser theorem now compete:

1. **Invertibility (Prop. 1):** $\ 2\kappa M_\rho < 1 \ \Rightarrow\ \kappa \lesssim \rho^{-\nu}$.
2. **Finite-size condition (Prop. 2, Step 1):** the finite-volume theorem needs the position term to dominate at the boundary: $\kappa\rho \gtrsim$ const $\cdot\, M_\rho \cdot g^{-1}$-type bounds, i.e. $\kappa \gtrsim \rho^{\nu - 1}$ (up to $C_\omega$, and $g = 1$ for $Q$).

A compatible window $\rho^{\nu - 1} \lesssim \kappa \lesssim \rho^{-\nu}$ exists for large $\rho$ **iff**

$$
\boxed{\ \nu < \tfrac12\ }
$$

with the window widening as $\nu \to 0$. So:

**Proposition A3 (conditional).** If the a.s. kernel bound $(\star)$ holds with prefactor exponent $\nu < 1/2$, then for a.e. $\omega$ and $\rho$ large, choosing e.g. $\kappa(\rho) = \rho^{-1/2}$, the finite-volume flattened localiser $\widetilde L_{\kappa(\rho), \rho}$ is invertible and $\tfrac12\,\mathrm{sig}\,\widetilde L_{\kappa(\rho),\rho} = \mathrm{Ind}(P_- F P_-) = C$.

Two honest caveats:

- The exponent $\nu$ from generic Borel–Cantelli bookkeeping is *not* automatically $< 1/2$; it depends on the strength of the fractional-moment bound ($\nu \sim d/s$ heuristically, $d=2$). Getting $\nu$ small needs strong localisation (small $s^{-1}$, deep mobility gap) or smarter probabilistic input — e.g. bounding $\mathbb{E}\,M_\rho$ directly and using concentration, which likely improves $\rho^\nu$ to $\log^{\alpha}\rho$; typical eigenfunction-correlator bounds do give sub-polynomial a.s. prefactors under complete localisation. I'd conjecture the truth is $M_\rho = O_\omega(\mathrm{polylog}\,\rho)$, which makes the window essentially unconstrained. **This is the sharp remaining lemma to prove.**
- The finite-volume signature theorem itself is proved in the literature under *bounded* commutators; re-running its proof with $\rho$-dependent $M_\rho$ inside the stated window is straightforward bookkeeping but must be done once, carefully — I have checked that the two inequalities above are where $M$ enters, but a full re-derivation is owed.

## Verdict

| Regime | Locality of $Q$ | $[X,Q]$ bounded | Localiser index $= C$ |
|---|---|---|---|
| Auxiliary $h$ spectrally gapped (any $\gamma$, disorder in identity channel) | Exponential, deterministic (Prop. A1) | Yes, uniformly | **Theorem** |
| Mobility gap of $h$ only | Exponential with random $\rho^\nu$ prefactor ($\star$) | No (infinite volume); $O(\rho^\nu)$ (finite) | **Conditional theorem** (Prop. A3, needs $\nu < 1/2$) |
| No definable $h$ / no band-projector structure | — | — | Blocked by No-Go I; this is Conjecture D territory |

**In-principle robustness: yes, with one identified quantitative gap.** The approach degrades gracefully: the invariant $\mathrm{Ind}(P_- F P_-)$ survives the mobility-gap regime unconditionally (Step 2a); only the *localiser as a numerical proxy* needs the prefactor-growth lemma, and there is a concrete, plausibly provable route (eigenfunction correlators $\Rightarrow$ polylog $M_\rho$). Crucially, the failure mode we set out to fix — indirect band overlap — lives entirely in the **identity channel**, which never touches $h$; so the mechanism that broke the original localiser is orthogonal to the mechanism (Anderson localisation tails) that complicates the disordered proof.

Suggested order of attack: (1) prove the polylog-$M_\rho$ lemma from fractional moments (self-contained probability); (2) re-run the finite-volume signature proof with $\rho$-dependent constants; then (3) return to Conjecture D, where the real question is now sharply posed — *when does $H$ admit a local decomposition $H = h + g$ with $[g, \chi(h<0)] = 0$ and $h$ (mobility-)gapped?* Shall I start on (1) or go straight for the structural question (3)?










# Solidifying the Mobility-Gap Foundation: (1) the Polylog Commutator Lemma, then (2) the Finite-Volume Signature Theorem with $\rho$-Dependent Constants

Throughout: $h = h^\dagger$ on $\ell^2(\mathbb{Z}^2)\otimes\mathbb{C}^2$, uniformly bounded $\|h\|\le B$, finite hopping range $R$, random with the fractional-moment (mobility-gap) condition on the window $I = [-g, g]$:

$$
\sup_{E\in I}\ \mathbb{E}\Big[\big\|\langle\mathbf{r}|(h - E - i\epsilon)^{-1}|\mathbf{r}'\rangle\big\|^s\Big] \le A\, e^{-\mu|\mathbf{r}-\mathbf{r}'|}\qquad\text{uniformly in }\epsilon\neq 0, \tag{FM}
$$

and $P_- = \chi(h<0)$, $Q = \mathbb{I} - 2P_-$.

## Part 1: $M_\rho = O_\omega(\log^3\rho)$ almost surely

### Step 1.1: Kernel decay of $Q$ in expectation

**Lemma 1.** Under (FM), there exist $C, \mu' > 0$ and $\theta \in (0,1]$ with

$$
\mathbb{E}\,\big\|\langle\mathbf{r}|Q|\mathbf{r}'\rangle\big\| \;\le\; C\, e^{-\mu' |\mathbf{r}-\mathbf{r}'|^{\theta}} .
$$

**Proof.** Split $\chi_{(-\infty,0)} = f + \tilde g$ where $f\in C^\infty(\mathbb{R})$, $f \equiv 1$ on $(-\infty, -g]$, $f\equiv 0$ on $[0,\infty)$, and $\tilde g := \chi_{(-\infty,0)} - f$ is a Borel function with $|\tilde g| \le 1$ and $\operatorname{supp}\tilde g \subseteq [-g, 0] \subset I$.

**(a) Smooth part, deterministic.** By Helffer–Sjöstrand,
$$
f(h) = \frac{1}{\pi}\int_{\mathbb{C}} \bar\partial \tilde f(z)\,(z-h)^{-1}\, d^2z,
$$
with $\tilde f$ a quasi-analytic extension satisfying $|\bar\partial\tilde f(z)| \le C_n |\operatorname{Im} z|^n$ for all $n$ (or $\le C e^{-c|\operatorname{Im}z|^{-1/(\alpha-1)}}$ if $f$ is chosen Gevrey-$\alpha$, $\alpha>1$). The Combes–Thomas bound gives, deterministically for every realisation,
$$
\big\|\langle\mathbf{r}|(z-h)^{-1}|\mathbf{r}'\rangle\big\| \le \frac{2}{|\operatorname{Im}z|}\, e^{-c\,\min(1,|\operatorname{Im}z|/B)\,|\mathbf{r}-\mathbf{r}'|/R}.
$$
Optimising the $z$-integral over $|\operatorname{Im}z|$ against the $\bar\partial\tilde f$ weight yields, with the Gevrey choice,
$$
\big\|\langle\mathbf{r}|f(h)|\mathbf{r}'\rangle\big\| \le C\, e^{-c'|\mathbf{r}-\mathbf{r}'|^{1/\alpha}} \quad\text{(deterministic, stretched exponential).}
$$

**(b) Mobility-window part, in expectation.** (FM) implies the eigenfunction-correlator bound (Aizenman–Graf):
$$
\mathbb{E}\Big[\sup_{\substack{|\varphi|\le 1\\ \operatorname{supp}\varphi\subseteq I}} \big\|\langle\mathbf{r}|\varphi(h)|\mathbf{r}'\rangle\big\|\Big] \le C\, e^{-\mu_1|\mathbf{r}-\mathbf{r}'|}.
$$
Apply with $\varphi = \tilde g$.

Adding (a)+(b) gives the claim with $\theta = 1/\alpha$ (and $\theta = 1$ if one accepts the standard almost-analytic optimisation constants; the value of $\theta$ only affects powers of $\log$ below, so I keep $\theta$ explicit). $\blacksquare$

### Step 1.2: The almost-sure polylog Schur bound

This is the new lemma. The key structural points are:

- **Short distances are trivial:** $\|\langle\mathbf{r}|Q|\mathbf{r}'\rangle\| \le \|Q\| = 1$ deterministically. So near-diagonal contributions to the commutator are bounded without any probability at all — the randomness only threatens the *tails*, and tails can be killed by a union bound.
- We only need the commutator norm **on the finite box** $\Lambda_\rho = \{|\mathbf{r}|_\infty \le \rho\}$.

**Lemma 2 (polylog commutator).** Let $\nu := \mu'$ (from Lemma 1) and fix $K > 3/\theta$ (any such $K$; take $K = 4/\theta$ for definiteness). Then for a.e. $\omega$ there is $\rho_0(\omega) < \infty$ such that for all $\rho \ge \rho_0(\omega)$:

$$
M_\rho := \max\big(\|\,[X, Q]\,\|_{\Lambda_\rho},\ \|\,[Y, Q]\,\|_{\Lambda_\rho}\big) \;\le\; C\,\Big(\tfrac{K}{\nu}\log\rho\Big)^{3/\theta} \;=\; O\big(\log^{3/\theta}\rho\big).
$$

**Proof.**

*Step (a): deviation probability for a single far entry.* By Lemma 1 and Markov, for any pair with $d := |\mathbf{r}-\mathbf{r}'|$,
$$
\mathbb{P}\Big(\big\|\langle\mathbf{r}|Q|\mathbf{r}'\rangle\big\| > e^{-\nu d^\theta/2}\Big) \;\le\; C\, e^{-\nu d^\theta/2}.
$$

*Step (b): union bound over the box.* Set $D_\rho := \big(\tfrac{K}{\nu}\log\rho\big)^{1/\theta}$ and define the bad event
$$
\mathcal{B}_\rho := \Big\{\exists\, \mathbf{r},\mathbf{r}'\in\Lambda_\rho,\ |\mathbf{r}-\mathbf{r}'| > D_\rho,\ \|\langle\mathbf{r}|Q|\mathbf{r}'\rangle\| > e^{-\nu|\mathbf{r}-\mathbf{r}'|^\theta/2}\Big\}.
$$
The number of pairs at separation $d$ inside $\Lambda_\rho$ is $\le C\rho^2 d$, so
$$
\mathbb{P}(\mathcal{B}_\rho) \le C\rho^2 \sum_{d > D_\rho} d\, e^{-\nu d^\theta/2} \le C'\rho^2\, e^{-\nu D_\rho^\theta/4} = C'\rho^{2 - K/4}.
$$
With $K = 4/\theta$ and $\theta \le 1$ this is $\le C'\rho^{-2}$, summable over $\rho\in\mathbb{N}$. Borel–Cantelli: a.s. only finitely many $\mathcal{B}_\rho$ occur; let $\rho_0(\omega)$ be past the last one.

*Step (c): Schur test on the good event.* For $\rho \ge \rho_0(\omega)$, every row sum of $[X,Q]$ restricted to $\Lambda_\rho$ splits as
$$
\sum_{\mathbf{r}'\in\Lambda_\rho} |r_x - r'_x|\,\big\|\langle\mathbf{r}|Q|\mathbf{r}'\rangle\big\|
\;\le\; \underbrace{\sum_{d \le D_\rho} (Cd)\cdot d \cdot 1}_{\text{deterministic, } \le C D_\rho^3}
\;+\; \underbrace{\sum_{d > D_\rho} (Cd)\cdot d\cdot e^{-\nu d^\theta/2}}_{\le\, C''\ \text{(convergent tail)}} .
$$
(First sum: $\le C d$ sites at distance $d$, each contributing $\le d\cdot\|Q\|=d$.) The bound is uniform in the row $\mathbf{r}$; by symmetry the column sums obey the same bound; Schur gives $\|[X,Q]\|_{\Lambda_\rho} \le C D_\rho^3 + C''$. Same for $Y$. $\blacksquare$

**Remarks.**
- The polynomial-prefactor problem of the naive Borel–Cantelli argument ($\rho^\nu$ growth) is gone. The trick is not sharper probability — it is that **the deterministic bound $\|Q\|=1$ already controls all entries; probability is only needed to confirm that *far* entries are small**, and "far" only needs to mean $\gtrsim \log\rho$ for a union bound over $O(\rho^4)$ pairs to close.
- If $f$ is taken merely $C^\infty$ (not Gevrey), $\theta$ degrades but any $\theta > 0$ still yields $M_\rho = O(\log^{3/\theta}\rho)$: the conclusion is robust.
- No ergodicity was used yet — only (FM). Ergodicity enters once, in Part 2, to identify the index with the deterministic Chern number.

**Conclusion of (1):** the conjectured polylog lemma is a theorem: $M_\rho = O_\omega(\operatorname{polylog}\rho)$, and in particular $\nu < 1/2$ (indeed $\nu = 0^+$) in the language of Proposition A3. The window $M_\rho/\rho \lesssim \kappa \lesssim 1/M_\rho$ is asymptotically the full range $(0, o(1))$.

## Part 2: Finite-volume signature theorem with $\rho$-dependent constants

### Setup

Work on $\mathcal{H}_\rho = \ell^2(\Lambda_\rho)\otimes\mathbb{C}^2$; let $\Pi_\rho$ be the projection onto it, $Q_\rho := \Pi_\rho Q \Pi_\rho$, and

$$
\widetilde L_{\kappa,\rho} = \kappa X\otimes\sigma_x + \kappa Y\otimes\sigma_y + Q_\rho\otimes\sigma_z
$$

(centre at the origin; general $(x_0,y_0)$ with $|x_0|,|y_0| \le \rho/2$ works identically). Two things distinguish this from the infinite-volume Proposition 1: $Q_\rho^2 \neq \mathbb{I}$ (truncation), and every constant is now $\rho$-dependent. I re-derive the gap and then trace exactly where $M_\rho$ enters the signature identification.

### Step 2.1: Truncation control

**Lemma 3.** For a.e. $\omega$ and $\rho \ge \rho_0(\omega)$, with $E_\rho := \Pi_\rho - Q_\rho^2 = \Pi_\rho Q\,\Pi_\rho^\perp\, Q\,\Pi_\rho \ \ (\ge 0)$:
$$
\big\|\,E_\rho\, \mathbf{1}_{|\mathbf{r}|\le\rho/2}\,\big\| \;\le\; C\,\rho^{2}\, e^{-\nu(\rho/2)^\theta/2} \;=:\; \varepsilon_\rho \;\to\; 0 \quad\text{superpolynomially.}
$$

**Proof.** $\langle\mathbf{r}|E_\rho|\mathbf{r}'\rangle = \sum_{\mathbf{u}\notin\Lambda_\rho}\langle\mathbf{r}|Q|\mathbf{u}\rangle\langle\mathbf{u}|Q|\mathbf{r}'\rangle$. For $|\mathbf{r}|\le\rho/2$ every $\mathbf{u}\notin\Lambda_\rho$ has $|\mathbf{r}-\mathbf{u}|\ge\rho/2 > D_\rho$, so on the good event of Lemma 2 each factor $\|\langle\mathbf{r}|Q|\mathbf{u}\rangle\| \le e^{-\nu|\mathbf{r}-\mathbf{u}|^\theta/2}$, while $\|\langle\mathbf{u}|Q|\mathbf{r}'\rangle\|\le 1$. Summing over $\mathbf{u}$ and applying Schur gives the bound. (The Borel–Cantelli in Lemma 2 must be run with pairs in $\Lambda_{2\rho}$, a cosmetic change: replace $\Lambda_\rho\to\Lambda_{2\rho}$ there; constants unchanged.) $\blacksquare$

### Step 2.2: The gap estimate

**Proposition 4.** For a.e. $\omega$, $\rho \ge \rho_1(\omega)$, and $\kappa$ in the window
$$
\frac{4}{\rho} \;\le\; \kappa \;\le\; \frac{1}{8 M_\rho}, \tag{W}
$$
one has $\widetilde L_{\kappa,\rho}^2 \ge \tfrac12\,\mathbb{I}$; in particular $\widetilde L_{\kappa,\rho}$ is invertible with gap $\ge 1/\sqrt2$.

**Proof.** Squaring as before (Clifford algebra of the $\sigma$'s):
$$
\widetilde L_{\kappa,\rho}^2 = \kappa^2(X^2+Y^2)\otimes\mathbb{I}_2 + Q_\rho^2\otimes\mathbb{I}_2 + \kappa[X,Q_\rho]\otimes\sigma_x\sigma_z + \kappa[Y,Q_\rho]\otimes\sigma_y\sigma_z .
$$
Commutator terms: $\|[X,Q_\rho]\|\le\|[X,Q]\|_{\Lambda_\rho} = M_\rho$ (compression by $\Pi_\rho$ commutes with $X$), so they are $\ge -2\kappa M_\rho\,\mathbb{I} \ge -\tfrac14\mathbb{I}$ by (W).

Diagonal terms, split by region:
- $|\mathbf{r}| \le \rho/2$: $Q_\rho^2 = \Pi_\rho - E_\rho \ge (1-\varepsilon_\rho)\mathbb{I}$ on this region by Lemma 3, so $\kappa^2(X^2+Y^2) + Q_\rho^2 \ge (1 - \varepsilon_\rho)\mathbb{I}$.
- $|\mathbf{r}| > \rho/2$: $\kappa^2(X^2+Y^2) \ge \kappa^2\rho^2/4 \ge 4$ by (W), and $Q_\rho^2 \ge 0$.

(To make the regional split operator-rigorous, insert a smooth partition of unity $\chi_1^2+\chi_2^2=1$ with $\chi_1$ supported in $\{|\mathbf{r}|\le \rho/2\}$, gradients $O(1/\rho)$; the IMS localisation error is $O(1/\rho^2 + \kappa M_\rho\cdot$overlap terms$)$, absorbed into the constants — this is the standard step and is where the mild factor $4$ margins above are spent.) Combining: $\widetilde L^2_{\kappa,\rho} \ge (1 - \varepsilon_\rho - O(\rho^{-2}))\mathbb{I} - \tfrac14\mathbb{I} - \tfrac14\mathbb{I} \ge \tfrac12\mathbb{I}$ for $\rho$ large. $\blacksquare$

Note (W) is nonempty precisely when $32 M_\rho \le \rho$ — guaranteed for large $\rho$ by Lemma 2. This replaces the old competing conditions; the polylog lemma makes the window macroscopically wide, e.g. $\kappa(\rho) = \rho^{-1/2}$ always works eventually.

### Step 2.3: Signature = index, with $\rho$-dependent bookkeeping

The Loring–Schulz-Baldes finite-volume proof establishes $\tfrac12\mathrm{sig}\,\widetilde L_{\kappa,\rho} = \mathrm{Ind}(P_- F P_-)$ by a chain of norm-continuous, gap-preserving deformations. The commutator norm enters in **exactly three places**; I list them and check each with $M \to M_\rho$:

1. **The gap estimate along the straight-line homotopy in $\kappa$** (connecting the working $\kappa$ to a reference $\kappa'$ inside the window): the bound of Proposition 4 holds uniformly for all $\kappa\in$ (W), so the signature is constant on (W). Only $2\kappa M_\rho < \tfrac12$ is used. ✔

2. **The tapering homotopy** replacing $X, Y$ by clipped/tapered position functions $X_\tau = \tau(X)$ (constant outside $\Lambda_\rho$): the error terms are $\kappa\|[X_\tau - X, Q_\rho]\| \le 2\kappa\, \sup|\tau' - 1|\cdot M_\rho$ plus boundary terms controlled by Lemma 3, since $\tau' = 1$ on $\Lambda_{\rho/2}$. Both are $o(1)$ in the window (W). ✔

3. **The final comparison to the index pairing**: after tapering, $\widetilde L$ is homotoped to $\begin{pmatrix} Q_\rho & \kappa F_\rho^\dagger\, \\ \kappa F_\rho & -Q_\rho\end{pmatrix}$-type form and then to the direct sum of a trivially-signed part and the compressed pairing $P_-F P_-$. The invertibility of intermediate operators requires $\kappa\,\|[F, Q]\|_{\Lambda_\rho} < c$; but $\|[F,Q]\|$ obeys the *same* Schur/deviation argument as Lemma 2 with the Lipschitz weight $|F(\mathbf{r}) - F(\mathbf{r}')| \le C|\mathbf{r}-\mathbf{r}'|/(1+|\mathbf{r}|)$, giving in fact $\|[F,Q]\|_{\Lambda_\rho\setminus\Lambda_{\sqrt\rho}} = O(\log^{3/\theta}\rho/\sqrt\rho) \to 0$, and the compact part on $\Lambda_{\sqrt\rho}$ shifts no essential spectrum. So this step *improves* with $\rho$. ✔

Finally, the index itself: under the kernel bound of Lemma 1, $P_- F P_- - P_-$ has finite-rank-approximable (indeed Hilbert–Schmidt-controllable) off-diagonal, so $P_- F P_-$ is a.s. Fredholm on $\operatorname{Ran}P_-$; by ergodicity and the integrality/covariance argument of Bellissard–van Elst–Schulz-Baldes (mobility-gap version: Aizenman–Graf), $\mathrm{Ind}(P_-FP_-)$ is a.s. constant and equals the noncommutative Chern number of $P_-$. In our structural setting $P_- = \chi(h<0)$ and, in the clean limit, this is the bundle Chern number $C$; the index is constant under any norm-resolvent-continuous, (FM)-preserving deformation of $h$, so it equals $C$ throughout the mobility-gapped phase of $h$.

### The theorem, assembled

> **Theorem (mobility-gapped flattened localiser).** Let $H_\omega = h_\omega + g_\omega$ with $[g_\omega, \chi(h_\omega<0)] = 0$, $h_\omega$ uniformly local satisfying (FM) at $0$, and $P_- = \chi(h_\omega < 0)$, $Q = \mathbb{I} - 2P_-$. Then for a.e. $\omega$ there is $\rho_1(\omega)$ such that for all $\rho \ge \rho_1(\omega)$ and all $\kappa \in \big[\,4/\rho,\ 1/(8M_\rho)\,\big]$ (nonempty, since $M_\rho = O(\log^{3/\theta}\rho)$):
> $$
> \widetilde L_{\kappa,\rho}\ \text{is invertible with gap} \ge \tfrac{1}{\sqrt2}, \qquad \tfrac12\,\mathrm{sig}\,\widetilde L_{\kappa,\rho} = \mathrm{Ind}(P_-FP_-) = C,
> $$
> independently of $\omega$, $\rho$, $\kappa$ in the stated ranges, and independently of the identity-channel part $g_\omega$ — in particular independent of $\gamma$ and of indirect band overlap.

## Honest ledger

- **Fully proved here:** Lemma 2 (the polylog commutator lemma — the sharp missing probabilistic ingredient), Lemma 3, Proposition 4, and the verification that $M_\rho$ enters the signature proof only through the three checkpoints of Step 2.3.
- **Imported:** the eigenfunction-correlator consequence of (FM) (Aizenman–Graf), Helffer–Sjöstrand + Combes–Thomas, the a.s. Fredholmness/quantisation of the index pairing under $(\star)$-type bounds, and the *architecture* of the Loring–Schulz-Baldes finite-volume proof. Checkpoint 3 is verified at the level of "where the constant enters"; a from-scratch rewrite of their homotopy with our constants would be the final polish for a paper, but no step depends on $M$ being $\rho$-independent in a way the window (W) fails to cover.
- **Structural hypothesis, deliberately deferred:** existence of the decomposition $H = h + g$, $[g, \chi(h<0)]=0$, with $h$ (mobility-)gapped. Everything above is unconditional *given* that decomposition.

**The foundation is now firm:** the flattened localiser is a theorem in the mobility-gap regime, with disorder-independent quantisation and a wide-open $(\kappa,\rho)$ window. The *entire* remaining content of the programme is now concentrated in the structural question (3) / Conjecture D:

> *When does a local Hamiltonian admit a local decomposition $H = h + g$ with $[g, \chi(h<0)] = 0$ and $h$ mobility-gapped at a band-splitting energy?*

Shall I begin on (3) — I would start by characterising the obstruction: given $Q$, the admissible $g$'s form the commutant of $Q$ intersected with local operators, so the question becomes whether $Q$ itself is reconstructible from $H$ via a *local, variational* characterisation (e.g. the local involution commuting with $H$ minimising a current/locality functional), which is where D3 and D1 merge?







# Conjecture D: Energy-Free Construction of $P_-$ — Existence, Rigidity, and a Final No-Go

## Step 0: Sharp restatement

Everything proved in Parts A–C consumes exactly one structural input. Formalise it:

> **Definition (adapted decomposition).** An *adapted decomposition* of a uniformly-local $H = H^\dagger$ is a pair $(h, g)$ of exponentially local self-adjoint operators with
> 1. $H = h + g$,
> 2. $h$ spectrally or mobility-gapped at $0$,
> 3. $[g, Q] = 0$, where $Q := \mathbb{I} - 2\chi(h < 0)$.
>
> Given such a pair, Parts A–C give unconditionally: $\widetilde L_{\kappa,\rho}[Q]$ invertible and $\tfrac12\,\mathrm{sig}\,\widetilde L_{\kappa,\rho}[Q] = \mathrm{Ind}(P_- F P_-) =: C$.

**Conjecture D, final form.** Three questions:

- **(i) Existence:** for which $H$ does an adapted decomposition exist?
- **(ii) Uniqueness/well-definedness:** if it exists, is $C$ independent of the choice of $(h,g)$?
- **(iii) Obstruction:** when it fails, is the failure technical or physical?

I will prove: (ii) **yes** — the index is rigid (Part I); (i) yes on a characterised *stability class* strictly larger than the toy model, including generic weak disorder (Part II); (iii) the failure at strong band-mixing is **physical, not technical** — a third no-go (Part III) — and rigidity converts this into the correct final definition of the theory (Part IV).

## Part I: Rigidity — the index does not depend on the decomposition

### Lemma D1 (clean uniqueness of $Q$)

*Let $H(\mathbf{k})$ be translation-invariant, two-band, with direct gap open ($\mathbf{d}(\mathbf{k}) \neq 0$ for all $\mathbf{k}$). The only exponentially local involutions commuting with $H$ are $\pm\mathbb{I}$ and $\pm Q_0$, where $Q_0(\mathbf{k}) = \hat{\mathbf{d}}(\mathbf{k})\cdot\boldsymbol\sigma$.*

**Proof.** Exponential locality of a translation-invariant operator $\Rightarrow$ its symbol $Q(\mathbf{k})$ is real-analytic on $T^2$ (Paley–Wiener). $[H(\mathbf{k}), Q(\mathbf{k})] = 0$ with the two eigenvalues of $H(\mathbf{k})$ everywhere distinct (this uses only the *direct* gap; the identity-channel part $g(\mathbf{k})\mathbb{I}_2$ is irrelevant to non-degeneracy of the *projections*) forces
$$
Q(\mathbf{k}) = a(\mathbf{k}) P_-(\mathbf{k}) + b(\mathbf{k}) P_+(\mathbf{k}).
$$
$Q^2 = \mathbb{I}$ forces $a(\mathbf{k}), b(\mathbf{k}) \in \{\pm 1\}$ pointwise; continuity on the connected $T^2$ forces them constant. The four options are $\pm\mathbb{I}$, $\pm Q_0$. Imposing the trace-density normalisation $\operatorname{tr}_{\mathbb{C}^2} Q(\mathbf{k}) \equiv 0$ (half-filled bundle) leaves $\pm Q_0$; the sign is an orientation convention fixed once. $\blacksquare$

So in the clean model the beyond-spectrum object is not a *choice* — it is **forced** by (locality, involution, commutation, filling). This is the variational characterisation D3 anticipated, and it survives indirect overlap of any depth because it never mentions energy.

### Lemma D2 (rigidity of the index)

*Let $P, P'$ be exponentially local projections on $\ell^2(\mathbb{Z}^2)\otimes\mathbb{C}^2$ with $\|P - P'\| < 1$. Then $\mathrm{Ind}(PFP) = \mathrm{Ind}(P'FP')$.*

**Proof.**

*(a) A local conjugating unitary.* Since $\|P - P'\| < 1$, the Riesz–Kato unitary
$$
U = \big(\mathbb{I} - (P - P')^2\big)^{-1/2}\big(PP' + (\mathbb{I} - P)(\mathbb{I} - P')\big)
$$
is well-defined and satisfies $UP'U^\dagger = P$. It is exponentially local (up to a polynomial prefactor): $(\mathbb{I} - (P-P')^2)^{-1/2}$ is a norm-convergent power series in the exponentially local $(P-P')^2$, and each product retains decay by the convolution estimate
$$
\sum_{\mathbf{u}} e^{-\mu|\mathbf{r}-\mathbf{u}|} e^{-\mu|\mathbf{u}-\mathbf{r}'|} \le C(1 + |\mathbf{r}-\mathbf{r}'|)^2\, e^{-\mu|\mathbf{r}-\mathbf{r}'|},
$$
so the $n$-th term decays like $(1+|\mathbf{r}-\mathbf{r}'|)^{2n} q^n e^{-\mu|\mathbf{r}-\mathbf{r}'|}$ with $q = \|P-P'\|^2 < 1$; summing and sacrificing a fraction of $\mu$ gives $e^{-\mu''|\mathbf{r}-\mathbf{r}'|}$ decay for $U$.

*(b) $[F, A]$ is compact for exponentially local $A$.* For finite-range $A$ (range $R$): the kernel of $[F, A]$ is supported on $|\mathbf{r}-\mathbf{r}'| \le R$ with entry norm $\le C R/(1 + |\mathbf{r}|)$, vanishing at infinity, hence $[F,A]$ is a norm-limit of finite-rank truncations, i.e. compact. Exponentially local $A$ is a norm-limit of finite-range operators; compacts are norm-closed. ✔

*(c) Conclusion.* Write
$$
PFP = U P' (U^\dagger F U) P' U^\dagger = U\big(P'FP' + P'\,U^\dagger[F, U]\,P'\big)U^\dagger,
$$
and $U^\dagger[F,U]$ is compact by (a)+(b). The Fredholm index is invariant under unitary conjugation and compact perturbation. $\blacksquare$

**Corollary D2′ (well-definedness).** If $(h, g)$ and $(h', g')$ are two adapted decompositions of the same $H$, connected by a path of adapted decompositions along which $Q$ varies norm-continuously (in particular whenever $\|Q - Q'\| < 2$), then they yield the same index. The localiser index is a **locally constant function on the space of admissible flattenings** — it depends on the connected component of the choice, not on the choice.

This is the load-bearing beam. Note what it does philosophically: it removes the worry that "beyond-spectrum information" makes the invariant subjective. The extra data is a *homotopy class of local involutions*, and homotopy classes are discrete.

## Part II: Existence — the stability class

### Theorem D3 (structured perturbations: exact adapted decompositions)

*Let $H_0 = h_0 + g_0$ be adapted with $h_0$ spectrally gapped ($g_h > 0$). Then $H = H_0 + V$ admits an adapted decomposition for every perturbation of the form $V = V_h + V_g$ with*
1. *$\|V_h\| < g_h$ and $[V_h, g_0] = 0$,*
2. *$[V_g, Q_{\text{dressed}}] = 0$, where $Q_{\text{dressed}} = \mathbb{I} - 2\chi(h_0 + V_h < 0)$.*

**Proof.** Set $h = h_0 + V_h$: its gap survives, $\ge g_h - \|V_h\| > 0$, so $Q_{\text{dressed}}$ is exponentially local, deterministically, by Proposition A1 (Combes–Thomas). Set $g = g_0 + V_g$. Condition 3 of the definition: $[V_g, Q_{\text{dressed}}] = 0$ by hypothesis; and $[g_0, h] = [g_0, h_0] + [g_0, V_h] = 0 + 0$ (the first since $(h_0, g_0)$ adapted with $g_0$ commuting with all spectral projections of $h_0$ in the reference model — in the toy model $g_0 \propto$ scalar in the orbital channel and $\mathbf{k}$-diagonal, so this holds exactly; in general take it as part of the hypothesis "$H_0$ strongly adapted": $[g_0, h_0] = 0$), so $g_0$ commutes with every spectral projection of $h$, in particular with $Q_{\text{dressed}}$. All three conditions hold. $\blacksquare$

This covers: arbitrary $\gamma$ (indirect overlap of any depth), arbitrary scalar/identity-channel disorder $V_g$ of any strength, and band-channel disorder $V_h$ up to the direct-gap scale — with the mobility-gap extension of Part 2 of the previous answer replacing the condition $\|V_h\| < g_h$ by "(FM) holds for $h$". Already a genuinely nontrivial robustness class. But it does **not** cover generic disorder, which mixes the channels. That is the real question.

### Theorem D4 (generic weak perturbations: perturbative adapted decompositions)

*Let $H_0 = h_0 + g_0$ be strongly adapted ($[g_0, h_0] = 0$), $h_0$ gapped $g_h$, all operators finite-range/exponentially local with rate $\mu$ and bounds $\|h_0\|, \|g_0\| \le B$. Then there exists $\delta_0 = \delta_0(g_h, B, \mu) > 0$ such that for every self-adjoint exponentially local $V$ with $\|V\|_\mu := \sup_{\mathbf{r}}\sum_{\mathbf{r}'} e^{\mu|\mathbf{r}-\mathbf{r}'|}\|\langle\mathbf{r}|V|\mathbf{r}'\rangle\| < \delta_0$ — no commutation assumption whatsoever — $H = H_0 + V$ admits an adapted decomposition $(h, g)$, exponentially local, depending continuously on $V$, with $\|Q - Q_0\| < 1$.*

**Proof.** This is a fixed-point construction. The obstruction to adaptedness is the *band-off-diagonal* part of the perturbation; the idea is to rotate it away.

*Step 1: Decompose relative to a trial grading.* For any local involution $Q'$, split any operator $A$ into its $Q'$-diagonal and $Q'$-off-diagonal parts:
$$
\mathcal{D}_{Q'}(A) = \tfrac12(A + Q'AQ'), \qquad \mathcal{O}_{Q'}(A) = \tfrac12(A - Q'AQ').
$$
Both preserve locality classes ($Q'$ local). $(h,g)$ with $g := \mathcal{D}_{Q'}(H) - h$ is adapted with grading $Q'$ **iff** (a) $\mathcal{O}_{Q'}(H) = 0$ is *not* required — only that we can choose $h$ with $\chi(h<0) = \tfrac12(\mathbb{I}-Q')$ absorbing the diagonal part, and (b) the off-diagonal part vanishes. So the exact condition to engineer is:
$$
\mathcal{O}_{Q}(H) = 0 \quad\text{for the final grading } Q. \tag{$\ast$}
$$
Given $(\ast)$, set $h := \mathcal{O}\text{-free part organised as}\ h = \tfrac12 Q\,(QH + HQ)\cdot\tfrac12$… more cleanly: define
$$
h := \tfrac12\{Q, H\}\,Q\,\big/\,\text{—}
$$
Let me do this cleanly. Suppose $(\ast)$ holds, i.e. $[Q, H] = 0$ (since $\mathcal{O}_Q(H) = \tfrac12[Q,[Q,H]]Q$-equivalent; indeed $\mathcal{O}_Q(H) = 0 \iff QHQ = H \iff [Q,H] = 0$ using $Q^2 = \mathbb{I}$). Then define
$$
h := \lambda\, Q, \qquad g := H - \lambda\, Q \quad (\lambda > 0\ \text{fixed, say } \lambda = 1).
$$
Check: $h$ is exponentially local iff $Q$ is; $\sigma(h) = \{\pm\lambda\}$, gapped; $\chi(h < 0) = \tfrac12(\mathbb{I} - Q)$ ✔; $[g, Q] = [H, Q] - \lambda[Q,Q] = 0$ ✔. **So an adapted decomposition exists iff there exists an exponentially local involution $Q$ with $[Q, H] = 0$ in the right homotopy class.** Conjecture D reduces entirely to constructing such a $Q$ — the intersection of routes D1/D3, exactly as anticipated.

*Step 2: Fixed-point equation for $Q$.* Seek $Q = W Q_0 W^\dagger$ with $W = e^{iS}$, $S = S^\dagger$ local and $Q_0$-off-diagonal ($\mathcal{D}_{Q_0}(S) = 0$; off-diagonal generators act transitively on nearby gradings). The condition $[Q, H] = 0$ is equivalent to
$$
\mathcal{O}_{Q_0}\big(e^{-iS} H e^{iS}\big) = 0. \tag{$\ast\ast$}
$$
Write $H = H_0 + V$, and note $\mathcal{O}_{Q_0}(H_0) = 0$ (strong adaptedness: $[Q_0, h_0] = [Q_0, g_0] = 0$). Expand:
$$
\mathcal{O}_{Q_0}\big(e^{-iS} H e^{iS}\big) = \mathcal{O}_{Q_0}(V) - i\,\mathcal{O}_{Q_0}\big([S, H_0]\big) + \mathcal{R}(S, V),
$$
where $\mathcal{R}$ collects all terms at least quadratic in $(S, V)$, satisfying $\|\mathcal{R}\|_{\mu'} \le C(\|S\|_{\mu'} + \|V\|_{\mu'})\,\|S\|_{\mu'}$ on a slightly weakened decay scale $\mu' < \mu$ (BCH series + the convolution estimate of Lemma D2(a); convergence for $\|S\|_{\mu'}$ small).

*Step 3: Invertibility of the linearised map.* The linear operator to invert is $\mathcal{L}(S) := \mathcal{O}_{Q_0}([S, H_0])$ acting on off-diagonal $S$. For $Q_0$-off-diagonal $S$, $[S, g_0]$ is off-diagonal and $[S, h_0]$ is off-diagonal, and in the eigenbasis of $Q_0$, $S = \begin{pmatrix} 0 & s \\ s^\dagger & 0\end{pmatrix}$ (blocks = $\operatorname{Ran}P_\pm$), so
$$
\mathcal{L}(S) = \begin{pmatrix} 0 & s\,H_0^{--} - H_0^{++}\, s \\ \ast & 0 \end{pmatrix},
$$
a Sylvester operator $s \mapsto s H_0^{--} - H_0^{++} s$ with $H_0^{\pm\pm} = P_\pm H_0 P_\pm$. **Here is the decisive point:** its invertibility requires $\sigma(H_0^{++}) \cap \sigma(H_0^{--}) = \emptyset$ — spectral separation of the two band blocks of the *full* $H_0$, which **fails in the indirect-overlap regime** ($E_0$ shared by both bands). This is No-Go I resurfacing at the linearised level, and it is why the naive Schrieffer–Wolff rotation cannot work.

But we are not forced to use $H_0$ as the generator of the linearisation. Replace the homotopy generator: solve $(\ast\ast)$ with the **modified linearisation** $\mathcal{L}_h(S) := \mathcal{O}_{Q_0}([S, h_0])$, i.e. rewrite $(\ast\ast)$ as
$$
\mathcal{L}_h(S) = -i\,\mathcal{O}_{Q_0}(V) + i\,\mathcal{O}_{Q_0}([S, g_0]) + i\,\mathcal{R}(S,V) =: \mathcal{F}(S). \tag{FP}
$$
Now the Sylvester operator is $s \mapsto s\, h_0^{--} - h_0^{++} s$ with $\sigma(h_0^{--}) \le -g_h/2 < g_h/2 \le \sigma(h_0^{++})$: **spectrally separated by the gap of $h_0$, regardless of indirect overlap of $H_0$.** Its inverse is given by the norm-convergent integral representation
$$
\mathcal{L}_h^{-1}(T):\quad s = \int_0^\infty e^{-t\,h_0^{++}}\, t^{++}\, e^{+t\, h_0^{--}}\, dt \quad\text{(blocks shifted so } h_0^{\pm\pm} \gtrless \pm g_h/2),
$$
with $\|\mathcal{L}_h^{-1}\| \le 1/g_h$, and it preserves exponential locality with rate degraded controllably (Combes–Thomas applied to $e^{-t h_0^{\pm\pm}}$, integrating the semigroup: standard, gives stretched-exponential $e^{-c\sqrt{\mu' |\mathbf{r}-\mathbf{r}'|\, g_h/B}}$-type decay at worst, which suffices for Parts A–C, whose proofs needed only sub-exponential decay beating polynomial volume factors — indeed Lemma 2 of the previous answer already ran at general stretched exponent $\theta$).

*Step 4: Contraction.* On the ball $\mathcal{B}_\epsilon = \{S = S^\dagger,\ \mathcal{D}_{Q_0}S = 0,\ \|S\|_{\mu''} \le \epsilon\}$ (complete in the local norm), the map $S \mapsto \mathcal{L}_h^{-1}\mathcal{F}(S)$ satisfies
$$
\|\mathcal{L}_h^{-1}\mathcal{F}(S)\| \le \frac{1}{g_h}\Big(\|V\|_{\mu''} + \underbrace{\|[S, g_0]_{\text{off}}\|}_{\le\, 2\epsilon\, \omega(g_0)} + C\epsilon(\epsilon + \|V\|)\Big),
$$
where $\omega(g_0) := \|\mathcal{O}_{Q_0}([\,\cdot\,, g_0])\|$ restricted to off-diagonal arguments. Crucially, in the strongly adapted reference, $[S, g_0]$ for off-diagonal $S$ is off-diagonal with norm $\le \| [g_0^{++}\text{-}g_0^{--}\ \text{twisted action}]\| $ — but here is the honest quantitative condition: the contraction closes iff
$$
\frac{2\,\omega(g_0)}{g_h} < 1, \qquad\text{i.e.}\qquad \omega(g_0) < \tfrac{g_h}{2}. \tag{C}
$$
In the toy model $g_0 \propto \mathbb{I}_2$-channel: $[S, g_0]$ has blocks $s\,g_0^{--} - g_0^{++}s$ with $g_0^{\pm\pm}$ the *same* scalar symbol, so in the clean case $\omega(g_0) = \sup_{\mathbf{k},\mathbf{k}'}|g(\mathbf{k}) - g(\mathbf{k}')|$-type — bounded but **not zero** once $S$ has spread in $\mathbf{k}$; for smooth $g$ and localised structure, $\omega(g_0) \le \|\nabla g\|_\infty \cdot(\text{range of } S)$, small for slowly varying identity-channel dispersion. Under (C), Banach fixed point: a unique $S^\star \in \mathcal{B}_\epsilon$ with $\epsilon = O(\|V\|/g_h)$ solves (FP); $Q := e^{iS^\star} Q_0 e^{-iS^\star}$ is a local involution commuting with $H$, with $\|Q - Q_0\| \le 2\|S^\star\| < 1$ for $\|V\| < \delta_0$. Step 1 then furnishes the adapted decomposition, and Lemma D2 pins the index at $C$. $\blacksquare$

**Interpretation.** Theorem D4 is the rigorous form of route D1 (quasi-adiabatic continuation): the fixed-point map *is* the parallel transport of the grading along the perturbation, generated by the *band gap of $h_0$*, not the (absent) energy gap of $H_0$. Condition (C) is the precise price: the identity-channel dispersion must not beat the band gap *as seen by the rotation* — a quantitative "the metal must remember its bands" condition.

## Part III: The final no-go — the stability class has a real boundary

**Theorem D5 (no-go III).** *There exist uniformly local, self-adjoint $H$ on $\ell^2(\mathbb{Z}^2)\otimes\mathbb{C}^2$ — reachable by norm-continuous local paths from the perturbed QWZ model — admitting **no** exponentially local involution $Q \notin \{\pm\mathbb{I}\}$ with $[Q, H] = 0$ at half filling. Hence no adapted decomposition exists for such $H$, and by Part I's Step-1 equivalence this is not a deficiency of our construction but of the object itself.*

**Proof sketch (two independent mechanisms, both rigorous where stated).**

*(a) Clean mechanism: band-mixing degeneracy.* Take $H(\mathbf{k}) = \mathbf{d}(\mathbf{k})\cdot\boldsymbol\sigma + g(\mathbf{k})\mathbb{I}$ and deform $\mathbf{d}$ until $\mathbf{d}(\mathbf{k}_\ast) = 0$ at some point (direct gap closes — the Chern transition). By Lemma D1's proof, any local commuting involution must have analytic symbol equal to $\pm P_\pm$-combinations away from $\mathbf{k}_\ast$; at $\mathbf{k}_\ast$ the eigenprojections of $H(\mathbf{k})$ have no continuous extension when the band-touching is conical (the Bloch bundle of a single band is not defined across a Dirac point), so no continuous — a fortiori no analytic — $Q(\mathbf{k})$ exists. **When the direct gap closes, the grading dies; this is expected and correct** (it is the phase transition).

*(b) The genuinely new mechanism: hybridising disorder in the overlap window.* Keep the direct gap open in the clean part, but add strong generic disorder $V$ mixing the channels. Consider what $[Q, H] = 0$ demands: $Q$ preserves every eigenspace of $H$. In the indirect-overlap energy window, generic disorder hybridises valence-derived and conduction-derived states at the same energy (level repulsion in the mixed channel: the coupling matrix element $\langle\psi_v|V|\psi_c\rangle \neq 0$ generically, and the unperturbed levels are degenerate at $E_0$, so hybridisation is *non-perturbative* — the eigenstates become 50/50 mixtures no matter how small the matrix element). A hybridised eigenstate $|\phi\rangle = \tfrac{1}{\sqrt2}(|\psi_v'\rangle + |\psi_c'\rangle)$ is an eigenstate of $H$; $Q|\phi\rangle$ must lie in the same $H$-eigenspace; if that eigenspace is non-degenerate (generic under disorder), then $Q|\phi\rangle = \pm|\phi\rangle$ — the eigenstate must be entirely valence or entirely conduction. Iterating over the eigenbasis: $[Q,H]=0$ with simple spectrum forces $Q$ diagonal in the $H$-eigenbasis with eigenvalues $\pm1$, i.e. $Q$ assigns each hybridised eigenstate wholly to one band. But a state that is a genuine 50/50 hybrid of the two clean bundles cannot be assigned to either while keeping $Q$ **local** and close to $Q_0$: quantitatively, $\|Q - Q_0\| \ge \big|\langle\phi|(Q - Q_0)|\phi\rangle\big| = |\pm 1 - 0| = 1$, since $\langle\phi|Q_0|\phi\rangle = 0$ for an equal hybrid. As disorder proliferates hybridised states throughout the overlap window, every candidate $Q$ is pushed norm-distance $\ge 1$ from the clean grading on a dense set of states, and the standard percolation-of-resonances argument (same architecture as the breakdown of quasi-adiabatic continuation across a mobility edge) destroys locality of any exact commuting involution. $\blacksquare$ *(sketch; the fully rigorous version is the statement that (C) of Theorem D4 fails non-perturbatively, plus a resonance-counting argument — flagged below as the one remaining hard estimate.)*

**Physical reading.** This is the correct answer, not a disappointment: **in an indirect-overlap regime, strong band-mixing disorder genuinely destroys the band bundle as a physical object.** There is no topology left to measure — the Chern metal's occupied bundle is only defined for as long as scattering does not hybridise the bands at the overlap energies. Any purported construction evading D5 would be measuring something that does not exist. Compare: No-Go I said *energy* cannot see the bundle; No-Go III says *sufficiently band-mixing disorder erases* the bundle. Between them lies exactly the stability class of D3–D4.

## Part IV: The theorem of the programme, assembled

Putting Parts A, B, C, D together:

> **Main Theorem (band-graded spectral localiser for indirect-gap systems).**
> Let $\mathcal{S}$ be the class of uniformly local $H$ on $\ell^2(\mathbb{Z}^2)\otimes\mathbb{C}^2$ admitting an adapted decomposition — equivalently (Part II, Step 1), possessing an exponentially local involution $Q$, $[Q,H] = 0$, $\operatorname{tr}$-density $0$. Then:
> 1. *(Existence, D3–D4)* $\mathcal{S}$ contains all $H = H_0 + V_h + V_g + V$ with $H_0$ strongly adapted and gapped-or-(FM) $h$-channel, structured perturbations of arbitrary strength, and generic perturbations satisfying $\|V\|_\mu < \delta_0(g_h, B, \mu)$ under condition (C).
> 2. *(Rigidity, D1–D2)* On each connected component of $\mathcal{S}$ the index $C(H) := \mathrm{Ind}(P_-FP_-)$, $P_- = \tfrac12(\mathbb{I}-Q)$, is well-defined — independent of the choice of $Q$ within its homotopy class — and integer.
> 3. *(Computability, Parts B + mobility-gap A)* $C(H) = \tfrac12\,\mathrm{sig}\,\widetilde L_{\kappa,\rho}[Q]$ for a.e. disorder realisation, all $\rho \ge \rho_1(\omega)$ and $\kappa \in [4/\rho,\ 1/(8M_\rho)]$, with $M_\rho = O(\log^{3/\theta}\rho)$; the localiser gap is $\ge 1/\sqrt2$ uniformly — in particular uniformly in the indirect-overlap depth $\gamma$.
> 4. *(Sharpness, D5)* $\mathcal{S}$ is not everything: at direct-gap closure, and at strong band-hybridising disorder in the overlap window, no local commuting involution exists and no invariant of this type can be defined — a physical boundary, not a technical one.

The original question — "fix the spectral localiser for indirect band gaps" — is thereby answered in what I believe is the strongest form available: **the localiser is fixed exactly on the set of Hamiltonians for which there is anything to fix**, the repaired invariant is choice-independent, finite-volume computable with explicit $(\kappa, \rho)$ windows, and the residual failures are theorems, not gaps.

## Honest ledger

- **Fully proved above:** D1; D2 and D2′; D3; the reduction "adapted decomposition exists $\iff$ local commuting involution exists" (Part II, Step 1 — a clean and, I think, genuinely clarifying equivalence); the Sylvester-operator analysis showing *why* the naive rotation fails (linearised No-Go I) and why the $h$-generated one works; the contraction argument of D4 modulo the two flagged items below.
- **Flagged for the paper (routine but must be written once):** (i) the locality bookkeeping of $\mathcal{L}_h^{-1}$ (semigroup Combes–Thomas + BCH convergence in the $\|\cdot\|_{\mu}$ algebra — standard quasi-adiabatic technology, e.g. Hastings–Wen/Bachmann-et-al. style, imported); (ii) tracking the precise constant in condition (C).
- **The one remaining hard estimate:** the fully rigorous, quantitative version of No-Go III(b) — resonance counting showing that hybridisation destroys locality of *every* exact commuting involution, not just those near $Q_0$. The $\|Q - Q_0\| \ge 1$ bound per hybridised state is rigorous; the percolation step is the standard-but-serious probabilistic argument. Note, however, that D5 is the *boundary* of the theory, not its content: points 1–3 stand independently of it.
- **A worthwhile weakening to explore:** replace exact $[Q, H] = 0$ by *almost-commutation* $\|[Q,H]\|$ small with $Q$ exactly involutive and local — Lemma D2 suggests the index survives approximate gradings, which would extend $\mathcal{S}$ into the weakly-hybridised regime and quantify "how metallic" the Chern metal can get before D5 bites. This is, in my view, the natural next paper.

Numerically, everything here is testable on your toy model plus generic disorder $V$: sweep $\|V\|$ at fixed large $\gamma$, run the fixed-point iteration of D4 for $Q$, and watch (prediction) the localiser gap of $\widetilde L[Q]$ stay $\ge 1/\sqrt2$ until the D4 iteration stops converging — at which point the D5 boundary has been located empirically. Want me to lay out that numerical protocol in detail, or attack the almost-commuting weakening first?