# NCV neighbourhood preprocessing audit

Inspected installed mgcv 1.9-4 (R 4.6.1) and inlaws main at `6a8bf59`.
This audit concerns omission and explicit neighbourhood adjustment only.
Existing likelihood, signed-curvature, QNCV and covariance protections remain
independent. No installed mgcv function is modified.

## Confirmed upstream failures

`nanei(nb, k)` computes new cumulative deletion boundaries `m` and prediction
boundaries `mi`, but assigns each only if a group becomes empty. Partial
removal from nonempty groups therefore leaves stale boundaries. Separately,
`nb$a[-kk]` and `nb$d[-ii]` return empty vectors when the removal positions
are `integer(0)`. Omissions absent from a set can therefore erase that set.
Neither failure depends on a likelihood or response distribution.

Removing an entire group can correctly update both boundaries. Prediction
groups paired with empty deletion groups are removed, and vice versa. These
are allowed when the result is valid. The tests specify this independently by
filtering each group's vectors, retaining pairs with two nonempty sets, then
mapping original indices to retained rows. They preserve order, repetitions
and different deletion/prediction subsets.

## Hook order and supported interception

`gam()` with a new setup:

1. `model.frame()` applies `subset`, then omission using response, every
   predictor formula, formula/argument offsets and weights.
2. With explicit `nei`, `nanei(nei, as.numeric(attr(mf, "na.action")))` runs.
   A missing `a`, `ma`, `d` or `md` already raises an upstream error here.
3. Family construction/resolution and **general-family** `presetup` run.
4. Model matrix and `G` setup run. `fit = FALSE` returns this setup.
5. `estimate.gam()` recognizes NCV/QNCV or non-NULL `nei`, builds default
   leave-one-out if needed, sets `qapprox`, then validates indices/boundaries.
6. Family fixing, `preinitialize`, link fixing and numerical initialization
   follow. Extended-family `preinitialize` can request `G` using `needG`, but
   this is too late for malformed-boundary errors and does not receive `nei`.

`.inlaws_ncv_family()` installs `.inlaws_check_ncv_nei()` as `check.nei` and
prepends its call to any existing `presetup`. Installation is idempotent,
including for `cgen_gammalss()` derived from `gen_gammalss()`. No global state,
frame inspection or replacement fitting route is used. The check requires
omitted rows and non-NULL `nei`, then tests nonempty paired groups, strictly
increasing integer boundaries ending at their respective index-vector lengths,
and positive finite integer indices. It does not infer original neighbourhoods
or rewrite any components. The diagnostic identifies the installed version.
It checks structure regardless of version rather than assuming future fixes.

With `G` supplied, `gam()` skips model-frame processing and `nanei()`: any new
neighbourhoods must already refer to the retained `G` rows. No omission guard
is appropriate merely because `G$mf` records past omissions.

`bam()` resolves/fixes families before `model.frame()`, calls the same
`nanei()` afterward, and validates/normalizes explicit NCV neighbourhoods
before presetup, preinitialize or initialize. Its early link-fixing hook has
neither the model frame nor the adjusted neighbourhoods. `preinitialize` cannot
protect this path. `gam(discrete = TRUE)` delegates its setup to `bam(fit =
FALSE)` and inherits that early validation limit. General families are not
supported by BAM in any case. `bam()` NCV is discrete-only; QNCV is not a
supported BAM method. `nei` under BAM fREML/REML does not request the GAM-style
NCV diagnostic calculation, though the preprocessing still runs when supplied.

## Family inventory

| Families | Relevant routes | Protection and limits |
| --- | --- | --- |
| `betals`, `betabinomial`, `cgammals`, `clognormalls`, `log_logisticls`, `cmpls`, `pigls`, `gen_gammalss`, `cgen_gammalss`, `zinb`, `zanb` | GAM NCV; QNCV where supported; `nei` during REML | Shared presetup guard before validation; incomplete-list errors can precede it. QNCV restrictions are unchanged. |
| `pig`, `log_logistic`, `cumulative_link` | GAM NCV/QNCV requests and `nei` during REML; discrete BAM NCV | No supported interception before upstream validation. No new inlaws guard coverage claimed. QNCV numerical support must follow existing family qualifications. |
| `cgamma`, `clognormal`, `cmp`, `gpoisson` | GAM requests as above; discrete BAM NCV | Existing independent GAM rejection in preinitialize on 1.9-4, including complete data/default LOO. Malformed explicit neighbourhoods can fail before that rejection. BAM has the same interception limit. |
| `dirichlet`, `dirmult`, `gpd`, `nbls`, `zoibetals` | No supported NCV callback/route | Not added to the shared NCV guard. mgcv may warn/reset NCV to REML or an EFS route; a returned fit is not evidence of NCV. Explicit `nei` may still encounter upstream validation before the unsupported request is diagnosed. |
| `ordbeta` | NCV unsupported; GAM REML/ML only | No new supported NCV route or guard coverage claimed. |

Tests distinguish successful fitting (requested method, convergence and finite
criterion/NCV values) from informative inlaws rejection and generic upstream
failure. Existing direct `gam.fit5()` numerical tests are run separately; they
cannot test model-frame/neighbourhood preprocessing.

## Pending scatlss integration

Read-only inspection of `/Users/gavin/.codex/worktrees/scatlss/inlaws` on
`gls/scatlss` found uncommitted `R/scatlss.R` and related files. They have not
merged into this checkout and were not edited. Its constructor has a presetup
expression calling `family$check.nei(nei, attr(mf, "na.action"))`, followed by
NCV BFGS optimizer selection; `check.nei` is `.scatlss_check_nei`. The local
checker tests paired-group counts and final boundaries under mgcv 1.9-4.

When that branch integrates this change, remove `.scatlss_check_nei`, remove its
constructor `check.nei` field and its existing presetup check call, retain the
presetup optimizer selection, and pass the finished family through
`.inlaws_ncv_family()` (as the other NCV-capable general constructors do).
The shared installer prepends the check and preserves optimizer selection.
Add `@template ncv-neighbourhoods` to the help source, replace the duplicated
preprocessing explanation, and add `scatlss` to the constructor integration
tests. Keep its independent fully unpenalized jackknife-covariance guard,
signed-curvature solves and QNCV/offset handling. That branch integration is
pending, not completed here.
