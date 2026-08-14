You: what is the existing literature on band topology. I know the Euler class for multi band topology, how does this relate? are there lessons from this existing area which can be bought into this localiser project. furthermore, are there results which already show/nullify the localiser results we are pursuing here?

Ask-Chat: Good moment to zoom out. I'll organise this as: (1) the map of existing band-topology theory and where this project sits, (2) the Euler class specifically, (3) transferable lessons, (4) prior art that directly overlaps — including work you *must* engage with because it's adjacent enough to threaten novelty — and (5) known results that constrain or could nullify the claims. One caveat up front: this field moves fast and my knowledge has a cutoff; for the 2024–26 literature you should enable web search and check recent arXiv listings (Cerjan, Loring, Schulz-Baldes, Slager, Bzdušek are the names to search).

## 1. The landscape

- **Energy-gap paradigm:** TKNN → Kitaev periodic table / Altland–Zirnbauer classes; K-theory of the Fermi projection (Bellissard–van Elst–Schulz-Baldes for disorder; Prodan–Schulz-Baldes monograph). Everything keyed to $\chi(H<E_F)$.
- **Real-space/local markers:** Bianco–Resta local Chern marker; Kitaev's local formula; the spectral localiser (Loring, Loring–Schulz-Baldes, finite-volume versions); Prodan's KPM-based invariants. All of these, in their original form, are energy-based — they inherit No-Go I.
- **Beyond the energy gap:** topological semimetals (Weyl/Dirac: invariants on spheres enclosing nodes), Fermi-surface topology (Haldane's AHE as Fermi-surface Berry phase), and "topological metals" broadly. Closest in spirit to your regime.
- **Multi-band / subbundle topology:** Wilson-loop spectra, Wannier obstruction theory, fragile topology (Po–Watanabe–Vishwanath; Cano et al.), topological quantum chemistry (Bradlyn et al.), and the Euler/Stiefel–Whitney class programme (Bzdušek, Bouhon, Slager, Ünal; Ahn–Park–Yang). This is where "topology of a *chosen set of bands* separated by direct gaps from the rest" is already the standard framing — in **k-space**.

Your project's slot in this map is quite specific: **real-space, disorder-compatible, finite-volume computation of subbundle (not Fermi-sea) invariants** — i.e. transplanting the multi-band k-space philosophy into the localiser/K-theory machinery. That intersection is much less populated than either parent area.

## 2. Euler class: how it relates

The Euler class $e \in H^2(T^2,\mathbb{Z})$ arises for **real** rank-2 bundles — in $C_2\mathcal{T}$ (or $\mathcal{PT}$)-symmetric systems where Bloch states can be chosen real. Key facts:

- It's the real-bundle analogue of the Chern class; nonzero Euler class obstructs Wannierisation of a two-band subspace even when all Chern numbers vanish.
- It's the paradigm of **fragile topology**: adding trivial bands to the subspace can kill it. Contrast: your Chern number $C(\mathcal{B})$ is stable (additive under enlarging $\mathcal{B}$, never destroyed by it).
- It governs non-Abelian braiding of band nodes between *different pairs* of bands (Wu–Soluyanov–Bzdušek; Bouhon et al.) — the "frame charge" picture where nodes carry quaternionic charges.

Relation to your construction: structurally it is exactly a "grading" invariant — it needs a direct-gap-separated two-band subspace, not a Fermi energy. So it lives on the same lattice of subbundle bipartitions your $\widetilde L[Q]$ explores. **But** your localiser cannot compute it as-is: the signature of $\widetilde L_\kappa$ computes the complex K-theory pairing (first Chern number). The Euler class is a KR-theory/real K-theory object; capturing it needs a localiser respecting the reality structure — there *is* existing work on spectral localisers in other symmetry classes (Loring–Schulz-Baldes cover the full tenfold way; Doll–Schulz-Baldes on skew localisers), so a "flattened Euler localiser" is a plausible follow-on project, but it's a different index pairing and fragility will fight the rigidity lemma (D2's homotopy argument is exactly what fragile invariants evade when you enlarge the space).

## 3. Lessons worth importing

1. **Gap conditions on band subsets are the standard currency** in the multi-band literature — your "direct-gap-separated grouping" is precisely their admissibility condition, so you can import their taxonomy (which groupings are protected, how invariants behave under subspace merging) wholesale.
2. **Wilson loops = the k-space shadow of your $Q$.** The Wilson-loop spectrum of a band group is computed from exactly $P_\mathcal{B}(\mathbf{k})$. Your localiser is the real-space, disorder-robust replacement for Wilson-loop winding. That's the cleanest one-line pitch of the project to that community.
3. **Fragile vs stable distinction:** be explicit that $\widetilde L[Q]$ captures stable (K-theoretic) invariants only. Rigidity (D2) is a feature for Chern numbers and a fundamental limitation for Euler-class physics.
4. **From semimetal theory:** invariants defined on subregions of the BZ (enclosing spheres) survive metallicity; Haldane's decomposition (quantised core + Fermi-surface correction) is the transport-side justification you're already leaning on — cite it as such.

## 4. Prior art that directly overlaps — engage or be scooped

This is the critical section. Several existing threads are close:

- **Cerjan–Loring and collaborators on localisers *without* spectral gaps.** They've shown the spectral localiser remains meaningful in gapless settings: the "local gap" of the localiser itself as an indicator, applications to Weyl/nodal semimetals and to photonic systems without bulk gaps, and the **quadratic pseudospectrum** (Cerjan–Loring–Vides). Their philosophy — "the localiser gap, not the Hamiltonian gap, is the protecting quantity" — partially anticipates yours. The distinction you must draw sharply: their constructions still put $H - E$ in the energy slot and study where its localiser gap happens to survive; yours replaces the slot's occupant with $Q$, which by No-Go I is not reachable from their framework in the indirect-overlap regime. That's a clean differentiation, but it must be argued explicitly against their papers.
- **Schulz-Baldes–Stoiber** on invariants of (disordered) semimetals and on localiser methods beyond gapped phases; also their work on approximate/mobility-gap regimes. Check whether their hypotheses cover indirect overlap — I believe they assume pseudo-gaps or vanishing DOS at $E_F$, which the Chern metal violates, but verify.
- **Band-flattening is standard in K-theory** ($H \to \mathrm{sgn}(H-E_F)$ is the first move in every proof). The nontrivial content of your work is flattening to a **non-spectral** grading. Someone skimming will say "flattening is old"; the answer is No-Go I — the flattening target is not a function of $H$. Make that the headline.
- **Quasi-adiabatic continuation** (Hastings–Wen, Bachmann–Michalakis–Nachtergaele–Sims): D1/D4 are essentially quasi-adiabatic transport of a *band* gap rather than an energy gap. The technology is imported; the application to gradings commuting with a metallic $H$ is, to my knowledge, new — but check Bachmann et al.'s later papers.
- **Almost-commuting matrices** (Loring's programme, Hastings–Loring): obstruction theory for "nearby commuting" involutions/projections is exactly the mathematical home of your No-Go III and of the proposed "almost-commuting $Q$" weakening. Results there (e.g. index obstructions to making almost-commuting matrices exactly commute) may *already contain* pieces of D5 — this could either strengthen the paper (cite and apply) or reveal the no-go is a known theorem in disguise.
- **Local markers for metals:** there's a strand (Bianco–Resta descendants; Marrazzo–Resta; more recent "Chern marker in metals" work) computing local Berry-curvature markers at fractional filling. These measure the Fermi-sea quantity (the non-quantised AHE piece), not your bundle invariant — again a differentiation to state, not a threat.

## 5. Could existing results nullify the claims?

Things that would sink or constrain it, as far as I know:

- **A theorem that any local involution commuting with a metallic local $H$ forces quantised transport** — would contradict the known non-quantisation of Chern-metal AHE. No such theorem exists; your framework avoids it precisely because $C$ is not claimed to be $\sigma_{xy}$. Consistent.
- **Non-existence results for local commuting structures:** Hastings-type results on the impossibility of localising Wannier functions for $C\neq 0$ don't obstruct you ($Q$ is an involution, not a frame — exponential locality of $P_-$ with nonzero Chern is fine; it's *Wannierisation* that's obstructed). Good — but say this, because referees will raise it.
- **The genuine risk is priority, not falsity:** the pieces (flattening, localiser-with-hypotheses-only, quasi-adiabatic band transport, mobility-gap pairings from Aizenman–Graf) all exist; the claim to novelty is the *assembly* plus No-Go I/III and the rigidity structure. I'm not aware of a published "band-projector localiser for indirect-overlap Chern metals", but given Cerjan–Loring's trajectory since ~2022 this is exactly where that programme is heading. **Do a serious literature sweep before investing in the paper** — search "spectral localizer metal", "spectral localizer band projector", "topological invariants indirect gap", "Chern metal invariant real space" on arXiv, with web search enabled since my information may be stale.

**Bottom line:** the Euler-class/multi-band literature validates your framing (subbundle topology with direct-gap admissibility is mainstream) but computes it only in k-space; the localiser literature has the real-space machinery but has stayed energy-anchored even in its gapless extensions. Nothing I know of nullifies the results; the closest work (Cerjan–Loring gapless localisers, Schulz-Baldes–Stoiber semimetals) must be engaged head-on to establish that "replace the energy slot with a constructed grading" is genuinely beyond their scope — which No-Go I argues it is.

You: can you do the web search right now please

Ask-Chat: I can't run a web search in this conversation — web search isn't currently enabled here, and it's something you control rather than me. To turn it on, look for the web search toggle in the chat input area, enable it, and then re-ask (e.g. "search for recent literature on spectral localisers for gapless/metallic systems").

Once it's on, here are the searches I'd suggest running, in priority order:

1. **`spectral localizer gapless metal arXiv`** — to catch the Cerjan–Loring line (local gap, quadratic pseudospectrum) and anything newer where they've moved toward metallic systems.
2. **`spectral localizer band projector flattened`** — the direct novelty check: has anyone already put a non-spectral grading in the energy slot.
3. **`topological invariant indirect band gap overlap Chern metal`** — for physics-side prior art on invariants surviving indirect overlap.
4. **`Schulz-Baldes Stoiber semimetal invariant`** — to pin down exactly what hypotheses their beyond-gap results assume (pseudo-gap / vanishing DOS vs genuine band overlap).
5. **`almost commuting matrices involution obstruction Loring Hastings`** — to see whether pieces of No-Go III already exist as theorems in the almost-commuting literature.
6. **`Euler class real space marker disorder`** — for the follow-on question of whether anyone has real-space Euler-class markers yet.

If you enable search and paste the question again, I'll go through the results and assess them against the specific claims in your document — particularly which items threaten priority on the $Q$-localiser assembly versus which are merely adjacent.
