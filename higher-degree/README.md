# Higher degree: the lifting problem in degrees 6, 7 and 8

This folder extends the degree ≤ 5 classification of
Kala–Yatsyna, *Even better sums of squares over quintic and cyclotomic fields*
(arXiv:2402.03850), to degrees 6, 7 and 8.

**Result.** No totally real field of degree 6, 7 or 8 admits a universal ℤ-form.
Together with the degree ≤ 5 theorem, the totally real fields of degree ≤ 8 with a
universal ℤ-form are exactly ℚ, ℚ(√5) and ℚ(ζ₇+ζ₇⁻¹).

The write-up with all proofs is [`notes/degrees6-7.pdf`](notes/degrees6-7.pdf)
(source [`notes/degrees6-7.tex`](notes/degrees6-7.tex)). The verification status
is listed in §9 of the notes (degree 8: §8.5–8.6).

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

The second paper (arXiv 0807.2099, hosted on hbs.edu) could not be downloaded from this
environment because both hosts are blocked, so it is not used here.

## Degree 8

**Case A (x₁ generates K): the rung split.** Let n = ⌈θ₁⌉ − 1, the lower rung, and
m = min(3ν/2, ν + 3/4). If τ(x₁ − n) < m, the square dichotomy makes x₁ − n = y² with
y ∈ O_K. Then one enumerates y instead of x₁: p₂(y) and p₄(y) are fixed, and |p₁(y)| is
small. The upper rung is symmetric. Otherwise both extreme conjugates lie a full unit
further out (θ₁ ≤ ⌊τ−m⌋+1, θ₈ ≥ ⌈τ+m⌉−1), and the tree pruning then collapses the tree.
`enum6.c` takes the arguments `X`, `YL` and `YU` for the three modes.

| step | count |
|---|---|
| leaves (all three modes, ν ≤ 3+2√2) | 6,012,361 (37 s on 4 cores; the plain run was days) |
| distinct polynomials / irreducible | 13,187 / 82 |
| exact tests in ℤ[x₁]: NS / BC / RG | 7 / 1 / 70 |
| maximal order (4 fields) | all fail sharpened rigidity |

Validation: on all 60 subtrees with S ≤ 24, every difference between the split run and
the plain run is either reducible or certified by `verify_split.gp` (99 + 5 and 46).
Degrees 6 and 7 were checked the same way.

**Case B (x₁ in a subfield F).**

* **Square forcing** (`sqforce8.gp`): every window element of F must be a square in K.
  * Of the 10 quartic F that pass NS/BC, 6 have ν(O_F) < Var(x₁). The other 4 force an
    octic field that still has a window non-square.
  * ℚ(√2) forces ℚ(ζ₁₆)⁺, then ℚ(ζ₃₂)⁺, which again has a window non-square.
  * ℚ(√5) forces nothing.
* **F = ℚ(√5)** (`rel8.gp`), with ν_{K/F} ≤ 4.21:
  * *Relative quartics* (`run_q4_par.sh`): 468,986 in the boxes. The exact relative
    downspread test (with the actual covering value for each target) removes 450,708.
    NSR/BCR removes 7,316, and 815 remain. Only 10 residues of c₁ mod 4 are needed, by
    √5 ↦ −√5.
  * *Intermediate quartic F′ = F(β)*: only d = 725 and ℚ(√2,√5) remain, and neither forces
    anything. The relative quadratics K/F′ (`enum_e2g8`, with an exact relative downspread
    prefilter) leave 52 + 97.
  * The 964 survivors give 390 fields. 386 fail relative sharpened rigidity (RGR2), 2 fail BCR,
    and 2 survive. Both fail condition (N) at α of trace 18 in the quartic subfield of
    discriminant 725 (minimal polynomial x⁴−9x³+27x²−31x+11):

| K | d_K | note |
|---|---|---|
| x⁸−2x⁷−12x⁶+26x⁵+17x⁴−36x³−5x²+11x−1 | 5⁴·29⁴ | Galois closure (D₄) of the quartic of discriminant 725 |
| x⁸−3x⁷−4x⁶+13x⁵+5x⁴−13x³−4x²+3x+1 | 5⁴·29²·1249 | non-Galois quadratic extension of that quartic |

`run_degree8.sh` reproduces everything. Case A takes about 1 minute. Case B takes about 2–3
hours on 4 cores, mostly the relative quartics and the maximal-order tests of 390 fields.

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

### Degree 8 checks (`validate8.sh`)

* **Case A exclusions re-checked exactly.** All 108,159 exclusions of irreducible
  polynomials made inside the split enumerator (NS, sharpened RIG, y-budget) were redone in
  exact arithmetic (`verify_split.gp`), with 0 disagreements. The other 484,381 excluded
  polynomials are reducible.
* **Condition (N) checked directly.** For all 392 field exclusions (Case A: 4 RG2;
  Case B: 386 RGR2, 2 BCR), `certN8.gp` decides (N) at the element the test names. (N)
  fails in every case, so each excluded octic field has a second, independent certificate.
* **Relative-quartic candidates confirmed by brute force.** `check_q4.gp` searches
  power-sum boxes and tests real-rootedness exactly (norm polynomial, Sturm). For all 16
  residues of c₁ it returns exactly the candidate set of `enum_q4`.
* **Split versus plain enumeration.** The two agree on all 115 subtrees of Case A
  (`compare_split.sh`), up to polynomials that are reducible or certified exactly: 144 found
  only by the plain run (137 reducible, 7 certified) and 46 only by the split run (all
  reducible).

## Files

```
run_degree6.sh          end-to-end reproduction, degree 6 (gcc, PARI/GP >= 2.15, python3+numpy)
run_degree7.sh          end-to-end reproduction, degree 7
run_degree8.sh          end-to-end reproduction, degree 8
run_sanity.sh           the same machinery in degrees 2-5 (answer key: Kala-Yatsyna)
validate.sh D           independent checks (exact re-check of C exclusions, enumerator)
src/enum6.c             Robinson/Rolle enumerator (compile with -DDEG=d); root filters and
                        tree pruning; mode "strong" adds the NS/RIG tests, "sharp" uses the
                        sharpened rigidity, "why"/"whysharp" print the reason for each exclusion
src/run_enum_par.sh     parallel driver for the enumerator
src/verify_strong.gp    exact re-check of the NS/RIG exclusions made in C
src/sharp.gp            sharpened rigidity (kappa >= 0): sqshape2, oshape2, fieldtest2
src/refined_bound.py    integrality-refined bound for nu per trace class
src/run_split_par.sh    degree-8 Case A: rung-split enumeration (modes X, YL, YU)
src/verify_split.gp     exact certificate for polynomials excluded only by the split
src/sqforce8.gp         Case B square forcing (subfields F, forced quadratic extensions)
src/rel8.gp             degree-8 Case B over Q(sqrt5): relative quartics, intermediate F',
                        relative quadratics over F' with exact relative downspread
src/run_q4_par.sh       parallel driver for the relative quartics
src/collect8.gp         maximal-order tests (fieldtest2) for the Case B survivors
src/witness_k1_d8.gp    (N) search for the degree-8 field of discriminant 5^4 29^4
src/certN8.gp           (N) decided at every degree-8 exclusion certificate
src/check_q4.gp         brute-force check of the relative-quartic candidates
src/compare_split.sh    rung split versus plain enumeration, differences certified
validate8.sh            the degree-8 checks
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

## Next: degree 9

The tensor floor, the square dichotomy and the rung split all still hold (tensor rank ≤ 8).
The new ingredient is cubic subfields with relative degree 3.
