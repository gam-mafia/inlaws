# Migration from mgcvUtils

The newly added families were extracted from the following local branch snapshots.
The original branches and worktrees are retained. Existing clognorm() and its
helper are not migrated. Smooth constructors are outside this package scope.

| Branch | Commit |
| --- | --- |
| `codex/betabinomial` | `44b0913633a544032ba54504a0664305a91d9430` |
| `codex/cgamma` | `0a30f66ac1791d6860784ccc40499173d85ebe6a` |
| `codex/clognormal` | `55ba09fd11fe1057e7a6fb6ff4d42942de0b36a2` |
| `codex/cmp-family` | `8d226ed6d4c08adb361c51fc1de8bbcb86786f88` |
| `codex/cumulative-link` | `531a17970f7c47658df77ddcab46d7a15738d2ad` |
| `codex/dirichlet-family` | `ebdf7e8730916ceef2a3e3af8ae0b3160d2ca659` |
| `codex/dirichlet-multinomial` | `32279966c43f803f99bdd66b0d09523e92e4aa2b` |
| `codex/generalized-poisson` | `7fd91e67c3304302a5af5e4872b08184c4090bc5` |
| `codex/ordbeta` | `727016d3ec0a2c181fd154173736edb8a76f5117` |
| `codex/zinb-zanb` | `1c8d98aac1246c36054d437fa423007a98d8a77e` |

Generalized Pareto code, help, and tests were copied from the uncommitted
working tree at `/Users/gavin/.codex/worktrees/gpd-family/mgcvUtils`
(branch `gls/gpd-family`, base `7f01f18`).

The `nbls` implementation comes from `codex/zinb-zanb`, which includes its
shared-derivative helper updates.

The package contains the families, their distribution callbacks, and the
family-specific `predict_ordbeta()` helper. Generic fitted-model simulation,
quantile, and family-fixing utilities remain outside inlaws. No dependency on
mgcvUtils or competing `gam` S3 methods is introduced.

Source license: MIT, copyright David L Miller (2020).

Tests and examples load inlaws and call distribution callbacks directly.
Checks specific to the original mgcvUtils wrappers (such as supplied-uniform,
extra-offset, and output-shaping behaviour) remain in the source branches;
family distribution checks are retained or compared with independent references.
The optional brms comparison remains in inst/validation, outside routine checks.

## Migration validation

Validated with R 4.6.1 and mgcv 1.9-4 on macOS (aarch64).

- All 21 revised family regression scripts and all documentation examples passed
  in `R CMD check --no-manual`.
- That run identified one namespace import note. After adding the explicit mgcv
  import, a check with `--no-manual --no-tests --no-examples` returned `Status: OK`.
  The import declaration was the only code-related change between those runs.
- Executable family code in all 17 family source files matches the source snapshots.
- The optional brms 2.23.0 Dirichlet likelihood and parameterization comparison passed.
- The local gratia `fix_family_rd()`, `fix_family_qf()`, and `fix_family_cdf()`
  implementations preserve all 34 supplied distribution callbacks and already
  provide the fallback distributions from the copied mgcvUtils utilities.
  Those utility functions belong in gratia; no additional functionality needed
  transferring for this migration. Gratia was inspected but not modified.
- No fitted-model simulation/quantile wrappers or family-fixing functions are
  included in inlaws. Tests call the family callbacks directly.

This validation does not establish support for every mgcv version or every
consumer of these callbacks. The original branches and worktrees remain intact.
No commit or push was made as part of the migration.
