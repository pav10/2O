# Higher degree: the lifting problem in degree 6

This folder extends the degree ≤ 5 classification of
Kala–Yatsyna, *Even better sums of squares over quintic and cyclotomic fields*
(arXiv:2402.03850), to degree 6.

**Result.** No totally real sextic field admits a universal ℤ-form. Together with
the degree ≤ 5 theorem, the totally real fields of degree ≤ 6 with a universal
ℤ-form are exactly ℚ, ℚ(√5) and ℚ(ζ₇+ζ₇⁻¹).

The write-up with all proofs is [`notes/degree6.pdf`](notes/degree6.pdf)
(source [`notes/degree6.tex`](notes/degree6.tex)). The verification status is
listed in §7 of the notes.

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

## What was computed

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

## Files

```
run_degree6.sh          end-to-end reproduction (gcc, PARI/GP >= 2.15, python3+numpy)
src/enum6.c             Robinson/Rolle enumerator (compile with -DDEG=d); root filters
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

## Next: degree 7

The theory holds verbatim for d = 7: tensor rank ≤ 6, and s/γ_s² ≥ 3/2 still
holds. A septic field has no proper subfields, so only Case A occurs, with
ν ≤ (11+√105)/4 ≈ 5.31. Compile `enum6.c` with `-DDEG=7` and run the order-level
filters.
