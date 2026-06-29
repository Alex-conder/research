# Response to Reviewers

**Manuscript**: Pirb-positive cells orchestrate post-ischemic neuroinflammation: a multi-dataset single-cell and spatial transcriptomics study  
**Authors**: [Author list]  
**Date**: 2026-06-16

---

We thank the reviewers for their constructive comments and suggestions. We have revised the manuscript accordingly. Below, we provide a point-by-point response to each comment.

## Reviewer #1

**Comment 1**: The study relies entirely on public datasets. Please clarify the rationale for dataset selection and discuss potential batch effects in more detail.

**Response**: We have added a new subsection in the Methods entitled "Dataset selection and batch-effect control" explaining the inclusion criteria (mouse MCAO/tMCAO models, post-ischemic time course, sc/sn/bulk/spatial modalities) and the ComBat-based batch correction strategy used for cross-dataset microglial integration. Limitations regarding residual batch effects are now explicitly discussed.

**Comment 2**: Can the authors exclude the possibility that Pirb+ astrocytes are contaminated by microglia doublets?

**Response**: We performed additional doublet filtering using Scrublet and recomputed the Pirb+ astrocyte DE analysis. The A1-like enrichment and top DE genes (Cd14, Fcer1g, C5ar1) remained significant, although we now emphasize that in situ validation (Pirb/GFAP/C3 RNAscope or IF) is required to confirm astrocytic identity. These results are included in Supplementary Table 2.

## Reviewer #2

**Comment 1**: The Visium data lack single-cell resolution. How confident can we be that Pirb+ spots reflect microglia rather than infiltrating monocytes?

**Response**: We agree that Visium cannot resolve individual cells. We therefore used microglial and inflammatory module scoring, correlated Pirb+ spots with known microglial markers, and interpreted the spatial data in conjunction with sc/snRNA-seq results. The revised Discussion clarifies this limitation and calls for future high-resolution spatial methods (e.g., Xenium, MERFISH) or combined immunofluorescence.

**Comment 2**: The functional validation roadmap is descriptive. Are there any preliminary experimental data?

**Response**: The current study is a multi-dataset computational and spatial analysis. The in vitro functional experiments are planned as a direct follow-up. We have made this distinction clearer in the Discussion and have included a detailed experimental design as Supplementary Methods.

## Reviewer #3

**Comment 1**: Please provide a more comprehensive reference list, especially regarding Pirb/LILRB2 in neuroinflammation and stroke.

**Response**: We have expanded the reference list to include key studies on PirB/LILRB2, myelin-mediated axon inhibition, reactive astrocyte biology, microglial heterogeneity in stroke, and spatial transcriptomics methods.

**Comment 2**: The Graphical Abstract is helpful but could better emphasize the peripheral-to-central transition.

**Response**: We have revised the Graphical Abstract to include an explicit arrow depicting Pirb+ peripheral monocyte/neutrophil infiltration into the ischemic brain and updated the legend accordingly.

---

We hope that our revisions adequately address the reviewers' concerns and that the manuscript is now suitable for publication in *Nature Communications*.
