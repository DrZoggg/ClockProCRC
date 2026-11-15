# Worked example: select HCC cohorts by analysis purpose

When comparing circadian-related cancer studies, record which dataset supports each analysis before comparing the reported results. This example turns the HCCcrds study's existing cohort descriptions into a filterable [five-dataset worksheet](data/hcc-cohort-roles.csv).

The source study is **Integrated-omics analysis defines subtypes of hepatocellular carcinoma based on circadian rhythm**, DOI [10.1016/j.joim.2025.06.003](https://doi.org/10.1016/j.joim.2025.06.003). It contributes circadian-related HCC subtyping and external subtype and survival comparisons. This is a cross-cancer evidence-comparison example, separate from the ClockProCRC scoring implementation.

| Dataset identifier as reported | Subtype comparison, Fig. 5A | Survival comparison, Fig. 5B |
| --- | --- | --- |
| ICGC_LIRI | Listed | Listed |
| GSE14520 | Listed | Listed |
| GSE54236 | Listed | Listed |
| GSE104580 | Listed | Not listed in this source extraction |
| GSE27150 | Not listed in this source extraction | Listed |

## Use the worksheet

1. **Choose the intended analysis.** For the reported survival comparison, filter the Fig. 5B column; for the reported subtype comparison, filter Fig. 5A. Each returns four listed datasets.
2. **Check overlap.** Three identifiers occur in both analysis sets. Their union contains five distinct dataset identifiers, rather than eight independent cohorts or publications.
3. **Keep reference experiments separate.** The article's circadian-reference experiments and single-cell framework construction serve different analytical roles from these external HCC comparisons. This worksheet extracts the selected Fig. 5 analysis sets; it is not the study's complete data inventory.

The worksheet preserves the original DOI, locator and source version in every row so it remains interpretable after download. Dataset identifiers are not a unique-participant count, and shared participants between source cohorts have not been adjudicated here.

## Source and version

The mapping reuses the already reviewed author-deposited, ahead-of-print text, **Results 3.5 and the Figure 5 legend**; final typeset figure-level equivalence has not been independently re-established for this worksheet. “Not listed” refers to this extraction, not to every possible use of that dataset elsewhere. The [source-linked author evidence page](https://drgezhang.com/papers/doi-10-1016-j-joim-2025-06-003.html#kf4) preserves the reviewed study context and original-source route.

For academic citation, cite the formal HCCcrds article DOI above. The worksheet is an author-maintained research aid, not another independent study or a new analysis result.
