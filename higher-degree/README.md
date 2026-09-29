# Higher degree: the lifting problem in degrees 6 and 7

This folder extends the degree ≤ 5 classification of
Kala–Yatsyna, *Even better sums of squares over quintic and cyclotomic fields*
(arXiv:2402.03850), to degrees 6 and 7.

**Result.** No totally real field of degree 6 or 7 admits a universal ℤ-form.
Together with the degree ≤ 5 theorem, the totally real fields of degree ≤ 7 with a
universal ℤ-form are exactly ℚ, ℚ(√5) and ℚ(ζ₇+ζ₇⁻¹).

The write-up with all proofs is [`notes/degrees6-7.pdf`](notes/degrees6-7.pdf)
(source [`notes/degrees6-7.tex`](notes/degrees6-7.tex)). The verification status
is listed in §8 of the notes.

## Why the degree ≤ 5 method does not extend, and what replaces it

In degree ≤ 5, Mordell–Ko turns a universal ℤ-form into "every element of 2O⁺ is a
sum of squares". That statement fails in rank 6 because of E₆. It also leads to a
house bound (1+√11) that is far too large to enumerate in degree 6.

Here the basic object is the **variance lattice** Λ_K = O_K/ℤ with the form Var.
Everything starts from condition (N) of [KY1, Prop. 3.1]: α = Σ M_ij ω_i ω_j with
M ⪰ 0 half-integral. From (N):

* **Tensor floor.** Tensor rank s ≤ 6 gives s/γ_s² ≥ 1, and ≥ 3/2 for s ≥ 2.
* **Budget–cost.** τ(α) = Tr(α)/d ≥ ν(K) = min Var on O_K∖ℤ, for every totally
  positive α ∉ ℤ.
* **Rigidity.** If τ(α) < (3/2)ν(K), then α = y² + cy + b with Var(y) ≤ τ(α).
  Equivalently, 4α − κ is a square in K for some κ ≡ 0, 3 (mod 4) in an explicit
  range, which is testable by factoring one degree-12 polynomial.
* **Relative versions** of the last two hold over every subfield F.
* **Downspread.** Budget–cost applied to x₁ − n and −x₁ − n′ gives
  ν(K) ≤ (5+√21)/2 for a sextic. Over a subfield F, it gives
  e·ν_{K/F} ≥ 2(ν_{K/F} − δ_F)², where δ_F is a "positive covering constant" of O_F.
* **E₆ threshold.** Even rank-6 lattices have det ≥ 3. Hence for
  Tr α < 3(3d_K)^{1/6} (≥ 29.47), (N) ⟺ 2α = Σ x_k² with Σ x_k ∈ 2O_K.
  This is used as an independent check.

## What was computed (degree 6)

| step | objects | outcome |
|---|---|---|
| A: shortest vector x₁ generates K | 2,088,667 real-rooted sextics with Var ≤ 4.79 → 40,764 pass the root filters → 15,720 irreducible | all 15,720 excluded inside ℤ[x₁] (NS 873, BC 5,062, RG 9,785); **no field computed** |
| B: x₁ in a subfield F | F ∈ {ℚ(√5), ℚ(√2), ℚ(ζ₇)⁺, ℚ(ζ₉)⁺}; relative enumeration | 62 order-level survivors → 18 fields |
| maximal-order tests | 18 fields | 12 fail relative rigidity, 1 fails rigidity |
| condition (N) | 5 fields | each fails (N) at an explicit α of trace ≤ 14 |

The five final certificates:

| K | d_K | witness α (min. poly) | Tr α |
|---|---|---|---|
| ℚ(√5, ζ₇⁺) | 300125 | x⁶−14x⁵+73x⁴−182x³+227x²−133x+29 | 14 |
| ℚ(√2, ζ₇⁺) | 1229312 | 2−√2 | 12 |
| ℚ(√2, ζ₉⁺) | 3359232 | 2−√2 | 12 |
| x⁶−9x⁴+24x²−17 | 7138368 | x³−6x²+9x−3 | 12 |
| x⁶−9x⁴+24x²−19 | 7978176 | x³−6x²+9x−3 | 12 |

## Degree 7

The same theory applies with ν(K) ≤ (11+√105)/4 ≈ 5.31. The tensor rank is ≤ 6, and
s/γ_s² ≥ 3/2 still holds. A septic field has no proper subfields, so only the case "x₁
generates K" occurs. To keep the output small, two order-level tests run inside the C
enumerator (mode `strong`):

* **NS:** c·x₁² + b·x₁ and x₁³ + a·x₁² + b·x₁ must not have smaller variance than x₁.
* **RIG:** rigidity for the downspread elements α±. This needs 4α − κ = z² in K. Since
  K = ℚ(α), the conjugates of z are ±√(4σⱼα − κ), and some sign pattern must give
  integral elementary symmetric functions.

| step | count |
|---|---|
| real-rooted septics with Var ≤ 5.31 (leaves) | 533,562,406 |
| excluded by downspread / quadratic / NS / RIG | 531,897,850 / 41,498 / 2,532 / 1,619,004 |
| remaining | 1,522 (23 irreducible) |
| exact tests in ℤ[x₁]: NS / BC / RG | 5 / 3 / 15 |
| **survivors** | **0**; no number field is computed |

`run_degree7.sh` reproduces this in about 1 minute on 4 cores, using the tree
pruning of §8 of the notes. The output is identical to the original ~30 minute run.

## Sharpenings (towards degree 8; §8 of the notes)

* **One universal form (via Oh, PAMS 128 (2000)).** For d ≤ 8, condition (N) holds for α
  iff α is represented over O_K by Φ_d = ½·Q of the even part of an n-universal lattice
  of minimal rank: ½(E₆⊥D₇) for d = 6, ½(E₈⊥D₇) for d = 7, ½(E₈⊥D₈) for d = 8. So K has a
  universal ℤ-form iff Φ_d is universal over K.
* **Sharpened rigidity.** In α = y² + cy + b, the constant κ = 4b − c² = Q(e)Q(r) − B(e,r)²
  is ≥ 0 by Cauchy–Schwarz. Hence:
  * every window element with τ(α) < min(3ν/2, ν + 3/4) is a **square in O_K**;
  * relatively, 4β − γ² is totally ≥ 0 in O_F.
* **Effects of sharpened rigidity:**
  * degrees 2 and 3 need no (N) at all;
  * in degree 6, 17 of the 18 Case-B fields die by rigidity, and only ℚ(√5, ζ₇⁺) needs (N);
  * in degree 6, F = ℚ(√2) is impossible in Case B, and F = ℚ(ζ₉)⁺ forces K = ℚ(ζ₃₆)⁺,
    which then fails rigidity.
* **Integrality-refined ν bound** (`src/refined_bound.py`): per trace class, ν ≤ 4.25
  (d = 6), 4.86 (d = 7), 5.50 (d = 8).
* **Tree pruning** (always on in filtered modes). Even centred moments satisfy
  M₄, M₆ ≥ aᵏ + bᵏ, and the middle roots lie in a window. The output is unchanged bit for
  bit, and the leaf count drops: degree 6 from 2.09M to 0.18M, degree 7 from 534M to 7.6M.
* **Degree 8 status:**
  * Case B, quartic F: 10 candidates. 7 are impossible, ℚ(√2,√5) is excluded, and
    ℚ(ζ₁₆)⁺ and ℚ(ζ₂₄)⁺ force K = ℚ(ζ₃₂)⁺ and ℚ(ζ₄₈)⁺, which both fail rigidity.
  * Case B, quadratic F: ℚ(√2) forces ℚ(ζ₃₂)⁺. **ℚ(√5) is the only open subcase**, needing
    an enumeration of relative quartics.
  * **Case A is still too expensive.** The subtree (a₁,a₂) = (0,−16) takes over 25 minutes,
    and the range goes down to a₂ = −20.

The second paper (arXiv 0807.2099, hosted on hbs.edu) could not be downloaded from this
environment because both hosts are blocked, so it is not used here.

## Checks

* `enum6.c` matches an exact-arithmetic re-implementation (`enum_check.gp`,
  Sturm sequences) on all 2,088,667 leaves in degree 6, and brute force in
  degrees 3–5.
* Run in degrees 2–5, the same pipeline recovers exactly the Kala–Yatsyna
  classification: ℚ(√5) and ℚ(ζ₇)⁺ survive, and everything else is killed. The
  (N) witnesses are 2−√2 in ℚ(√2) and 2−ζ₉−ζ₉⁻¹ in ℚ(ζ₉)⁺. Degree 4 leaves three
  fields, all failing (N).
* Every (N)-witness is confirmed twice: by the M-search (`condN.gp`) and by the
  sum-of-squares search below the E₆ threshold (`sos2.gp`).
* Every exclusion made by the in-enumerator tests (NS/RIG) on an irreducible
  polynomial was re-checked in exact arithmetic (`verify_strong.gp`). That is
  15,714 polynomials in degree 6 and 857,103 in degree 7, with no disagreements.
  In degree 7 the enumerator also agrees leaf by leaf with `enum_check.gp` on the
  48 subtrees with the fewest leaves (2,267,434 leaves).
  `validate.sh D` reruns these checks.

## Files

```
run_degree6.sh          end-to-end reproduction, degree 6 (gcc, PARI/GP >= 2.15, python3+numpy)
run_degree7.sh          end-to-end reproduction, degree 7
run_sanity.sh           the same machinery in degrees 2-5 (answer key: Kala-Yatsyna)
validate.sh D           independent checks (exact re-check of C exclusions, enumerator)
src/enum6.c             Robinson/Rolle enumerator (compile with -DDEG=d); root filters and
                        tree pruning; mode "strong" adds the NS/RIG tests, "sharp" uses the
                        sharpened rigidity, "why"/"whysharp" print the reason for each exclusion
src/run_enum_par.sh     parallel driver for the enumerator
src/verify_strong.gp    exact re-check of the NS/RIG exclusions made in C
src/sharp.gp            sharpened rigidity (kappa >= 0): sqshape2, oshape2, fieldtest2
src/refined_bound.py    integrality-refined bound for nu per trace class
src/run_enum.sh         runs the enumerator over all (a1,a2)
src/enum_check.gp       exact re-implementation of the enumerator, for validation
src/lift6.gp            variance lattice of an order: nu, totally positive elements of small
                        trace, rigidity search, rigidity-via-squares test (sqshape)
src/order_filter.gp     order-level tests NS / BC / RG for Case A
src/deltaF.py           certified upper bounds for the positive covering constants delta_F
src/rel6.gp             Case B: relative enumeration over F and order-level tests in O_F[y]
src/field6.gp           maximal-order tests: BC, RG and relative BC/RG over all subfields
src/condN.gp            exact decision procedure for condition (N)
src/sos2.gp             independent check: sums of squares with parity + E6 threshold
src/witness.gp          search for an (N)-witness by increasing trace
src/pipeline.gp         drivers used by run_degree6.sh
results/                all outputs (status of every candidate, survivors, witnesses)
notes/                  LaTeX write-up
```

## Next: degree 8

See "Sharpenings" above. The only open Case B subcase is F = ℚ(√5), which needs relative
quartics with ν_{K/F} ≤ 4.21. Case A needs further pruning before it is practical.
