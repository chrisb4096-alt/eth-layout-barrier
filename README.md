# Lower bounds for time-dependent layout measurements

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23845341.svg)](https://doi.org/10.5281/zenodo.23845341)

Chris Blaschke (Ax0s), October 31, 2026. [Read the paper (PDF).](paper/layout-barrier.pdf)

This paper proves lower bounds for the accounting exponent in the layout interface of OpenAI's *The Exponential Time Hypothesis is False*. The bounds apply to bounded-incidence, balanced-cut expanding scopes and allow the measurement spaces to change at every evaluation point.

## What is new

The earlier bounds of [Tulsyan](https://github.com/xennygrimmato/eth-hamiltonicity/blob/a4f924c1419cbc43564c3792925d8bb182ca932c/research/permanent/parity-barrier.md) and [Olivier](https://gist.github.com/antoine-olivier/d269b3b449d3a2cdb021dc0ee5d0fc82) apply to the fixed-$M_0$ template, not to every layout. Tulsyan distinguishes this template from the general layout interface and leaves a lower bound for general time-dependent measurements open, noting that a proof would have to couple corrections that differ at each evaluation point. We couple them through one determinant: rank drops at different evaluation points become roots of a single polynomial of degree at most $mR$. This gives:

- an inequality $m(L-R)\le\sum_t c_t$ for every layout, with no expansion or integrality assumption (Theorem 1);
- the bounds below on expanding scopes (Theorem 3);
- in Tulsyan's auxiliary coordinates, on the full label space and for every time-dependent correction, $R\ge c\min\{Q,\sqrt{Q\ln(1/\varepsilon)}\}$ when $\varepsilon\le K/Q$ (Theorem 4);
- a second argument through interpolation at "partner" times, which is tight on a family of binary coset instances whose full-set row and column group counts have zero imbalance (Proposition 5, Corollary 6).

## Main bound

Suppose each variable occurs in at most a fixed number $b$ of scopes, and every set $W$ with $n/4\le |W|\le n/2$ is crossed by at least $h|W|$ factors, for fixed $h>0$. For any layout, with its parameters as defined in the paper, let

$$
F=\frac{R}{Q}\log_2\frac{2Q}{R}+\tau\log_2(2D+1)+\frac1r.
$$

There is a constant $c>0$ depending only on $b,h$ such that

$$
F\ge c\frac{\log_2(2Q)}{\sqrt Q}.
$$

If $Q\log_2 r\le K\log_2 n$ also holds for fixed $K$, then for sufficiently large $n$,

$$
F\ge c_{b,h,K}\frac{(\log_2\log_2 n)^{3/2}}{\sqrt{\log_2 n}}.
$$

The expansion, bounded-incidence and capacity assumptions are not part of the layout definition. The result bounds the upper-bound accounting expression $F$ used in OpenAI's analysis for total product degree $O(n)$, not the actual running time of every algorithm. It does not prove the capacity-free bound or match the paper's full logarithmic order.

## Verify

From the repository root, with `uv` installed:

```sh
uv run --with sympy python verify/check_timedep.py
uv run --with sympy python verify/partner_identity.py
```

The first command checks finite layout, determinant and probability certificates and writes `verify/certificates.json`. The second checks a projection identity and block-rank identity. These checks do not prove the general theorems.

To compile the paper with Tectonic:

```sh
cd paper
tectonic -X compile layout-barrier.tex
```

## Citation

```bibtex
@misc{blaschke2026layout,
  author = {Chris Blaschke},
  title = {Lower bounds for time-dependent layout measurements},
  year = {2026},
  month = {October},
  note = {Ax0s. First posted October 31, 2026},
  doi = {10.5281/zenodo.23845341},
  url = {https://doi.org/10.5281/zenodo.23845341}
}
```

## License

Paper text: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Verification code: [MIT](https://opensource.org/license/mit).

First posted October 31, 2026.
