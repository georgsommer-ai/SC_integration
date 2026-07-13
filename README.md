Author: Georg Sommer

# Disentangling Cell Identity and Interferon Response in Human PBMCs

**Integrated UMAP of control and interferon-stimulated PBMCs.**  
CCA-based integration reduces treatment-driven separation while preserving biologically meaningful cell-type structure. Control and stimulated cells occupy comparable cell-type neighborhoods, enabling cell-type-specific comparison of the interferon response.

## Project overview

Interferon stimulation induces broad transcriptional changes across peripheral blood mononuclear cells (PBMCs). These changes can become so dominant that cells cluster according to treatment condition rather than biological cell identity.

This project uses a Seurat-based single-cell RNA-sequencing workflow to separate treatment-associated expression changes from stable PBMC cell identities and to characterize the interferon response of various PBMCs.

## Research question

**Can integration recover biologically meaningful PBMC cell identities despite a strong interferon-induced transcriptional shift?**
**What patterns do we see in the response of PBMCs to interferon stimulation?(Ctrl vs. Stim in the variables and plots)**

## Summary

CCA-based integration aligned control and interferon-stimulated PBMCs by cell identity, after which differential-expression and visualization analyses revealed a broad interferon-response program.

## Aims

This project aims to:

- determine how strongly interferon stimulation changes the transcriptional structure of the PBMC dataset;
- distinguish stable cell-type identity from treatment-associated variation;
- compare selected interferon-response genes across PBMC cell types.

## Dataset

The analysis uses the `ifnb` dataset distributed through the `SeuratData` package.

The dataset contains approximately:

- **14,000 human PBMCs**
- **two experimental conditions**
  - `CTRL`: untreated control cells
  - `STIM`: interferon-stimulated cells
- multiple annotated immune-cell populations

## Analysis strategy

The workflow first creates an unintegrated representation of the data. This representation is used as a diagnostic step to determine whether treatment-associated variation dominates the analysis.

The treatment groups are then integrated using canonical correlation analysis. Clustering and UMAP are repeated in the integrated space to obtain a representation that better reflects shared cell identities.

The interferon responses of different PBMCs subpopulations are compared with Differential expression.

## Differential expression results

**Are the top interferon response markers different in CD14 and CD16 monocytes?**
There is minimal difference in interferon response markers in CD14 and CD16 monocytes.

**Does activation of B cells influence B cell interferon response?**
Activation of B cells has minimal influence on B cell interferon response

![Figure 1: B cell state has little to no influence on interferon response](docs/figures/B_response.png)

<br>
<br>
<br>

**What influence does T cell type and state have on the interferon response?**
T cell type and state have little to no influence on interferon response

![Figure 2: CD4 T cell state has little to no influence on interferon response](docs/figures/CD4_response.png)

<br>
<br>
<br>

## Key takeaways

- Interferon stimulation produces a broad transcriptional response across human PBMC populations.
- Unintegrated clustering can be dominated by experimental condition.
- CCA-based integration improves alignment of corresponding cell types across conditions.
- Stable lineage markers and treatment-responsive genes can be distinguished after integration.
- `ISG15` and `CXCL10` are prominent examples of the general interferon-response.
- PBMCs show a measurable condition-specific response that can be investigated using differential expression.

## Limitation

- Strong biological conclusions require appropriate replication, quality control, statistical testing, and independent validation.
