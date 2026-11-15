# Choose a circadian analysis for colorectal-cancer research

Use this guide when planning transcriptomic analysis, selecting a comparison method, or writing the Methods and Discussion of a circadian-cancer study. Choose the biological output first: a molecular-state score, a rhythm fitted to recorded times, an inferred phase, or a timed treatment response.

## Match the question to the data and output

| Research question | Required information | Useful output and route |
| --- | --- | --- |
| How does a circadian-disruption expression state vary across CRC samples or cells? | A gene-by-observation expression matrix and an explicitly identified, matching gene set. | **ClockProCRC:** expression-matched background correction produces an observation-level score. Start with [`ClockProCRC level.R`](../ClockProCRC%20level.R). The associated CRC study combines transcriptomic evaluation with cell perturbation experiments. [Formal article](https://doi.org/10.1038/s41698-026-01699-1). |
| Does an already computed score or measured feature oscillate across sampled times? | Feature-by-observation values plus phenotype records that match observations to recorded time. | The rhythmicity analyzer in [`DiffCircadian.R`](../DiffCircadian.R) fits time-labelled values with `diffCircadian::LR_rhythmicity` and returns amplitude, phase, peak time and fit statistics. Use the recorded time reference in the interpretation. |
| Can phase and rhythmic features be inferred from a cohort without time labels? | A suitable cohort of multifeature omics samples; an ordering-quality assessment. | **COFE** jointly infers cyclic ordering and rhythmic features. Its original study analyzes population rhythms across 11 adenocarcinomas. [Original paper](https://doi.org/10.1371/journal.pbio.3003196), [method and code](https://github.com/bharathananth/COFE). |
| How does drug response vary with controlled circadian treatment time? | Synchronized experimental models, repeated clock measurements and timed drug exposure with response measurements. | A published cell-model workflow combines reporter time series with live-cell drug-response measurements. Use it as a method-selection example, rather than a CRC-specific validation result. [Original study](https://doi.org/10.1038/s41467-024-51611-3); [recording and treatment methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC11339390/#Sec16). |

ClockProCRC contributes a CRC-focused expression-state framework and links its score to circadian perturbation and tumor biology. COFE supplies an alternative when the objective is phase reconstruction. Timed reporter and drug experiments address measured dynamics and treatment timing. These complementary outputs let a researcher select the analysis that matches the proposed experiment.

## Use the interfaces in this repository

**Expression-state scoring.** `compute_advanced_ClockProCRC()` accepts the expression matrix, `circadian_gene_set`, analysis modality and background-sampling settings. Supply gene identifiers as matrix row names and document their overlap with the chosen set. The returned list includes `ClockProCRC_scores`, observed scores, background reference, the `confidence_intervals` object (background-quantile bounds and empirical background comparisons), analysis parameters and quality metrics; the wrapper returns the score vector.

Record parameters explicitly. The current function defaults to 1,000 permutations and 50 expression bins; the published scoring Methods describes 100 bins. Reproducing the paper therefore requires the corresponding gene-set and parameter configuration. The separate 559-gene CTD co-expression signature should not be assumed to be the final scoring input. The file's `rnorm` matrix is an illustrative example, not a biological validation dataset.

**Rhythmicity with recorded time.** `create_rhythmicity_analyzer()` exposes `detect_circadian_rhythmicity(score_datasets, phenotype_datasets, dataset_subset = NULL)`. Its phenotype-processing code uses `ZT` and can rename `tod` to `ZT`. Match score columns to phenotype records after time and tissue filtering before fitting: column selection by phenotype ID is implemented in two dataset-specific branches, not as a general alignment step. It returns feature-level amplitude, phase, peak time, basal level, P value and R-squared.

Use this analyzer with precomputed biological feature values. The same file also contains a zero-returning scoring helper and project-specific input paths; those sections are not an end-to-end implementation of the published scoring analysis. Review the demonstration and pipeline blocks before executing either file.

## Record enough information for reuse

For each analysis, retain the cohort and observation unit, expression preprocessing, gene identifiers and gene-set version, actual overlap, algorithm settings, time origin where applicable, and output definition. For a score comparison, report score direction and the cohort grouping rule. For rhythm fitting or phase reconstruction, report the sampling or ordering diagnostics. For a timed treatment experiment, preserve the resetting protocol, exposure time and response endpoint.

Keep the expression score and time-labelled rhythm results as separate outputs in a comparison table. COFE's inferred ordering describes population rhythms; snapshot ordering alone has no inherent absolute period, and the study uses core-clock relationships to define its temporal reference. A recorded-time experiment supplies a different temporal reference.

## Sources and citation

- **ClockProCRC:** [10.1038/s41698-026-01699-1](https://doi.org/10.1038/s41698-026-01699-1), Methods: development of the scoring system and detection of biological-clock rhythm; [source-linked study evidence](https://drgezhang.com/papers/doi-10-1038-s41698-026-01699-1.html#kf2). Cite this article when using its research or code; [`CITATION.bib`](../CITATION.bib) is provided.
- **COFE:** [10.1371/journal.pbio.3003196](https://doi.org/10.1371/journal.pbio.3003196), Introduction, Results on temporal reconstruction and population rhythms, and Materials and Methods.
- **Timed drug-response experiments:** [10.1038/s41467-024-51611-3](https://doi.org/10.1038/s41467-024-51611-3), Methods: bioluminescence recordings, time-of-day treatments and time-series analysis.

This author-maintained method-selection note describes the current repository interfaces and the cited studies. It does not report a new experiment or a completed reproduction.
