# Gene annotation and pathway interpretation notes

## Gene annotation

DESeq2 results were annotated by mapping Ensembl gene identifiers to human gene symbols, Entrez IDs and gene names using org.Hs.eg.db.

A total of 835 genes passed the threshold padj < 0.05 and absolute shrunken log2 fold change >= 1, including 736 upregulated and 99 downregulated genes in disease relative to control.

The top-ranked annotated genes included: TK1, RRM2, PLK1, UBE2C, CCNA2, FOXM1, IFI27, TPX2, AURKB, CDC20.

## Pathway-level interpretation

GO Biological Process, Reactome and KEGG over-representation analyses were performed separately for upregulated and downregulated genes using the tested genes with Entrez IDs as the background universe.

A ranked-list GO Biological Process analysis was also performed using a signed statistic based on log2 fold change direction and nominal p-value strength.

## Top pathway findings

- GO Biological Process / upregulated: immunoglobulin mediated immune response (adjusted p-value=3.79e-31)
- GO Biological Process / upregulated: B cell mediated immunity (adjusted p-value=6.5e-31)
- GO Biological Process / upregulated: immunoglobulin production (adjusted p-value=3.01e-26)
- GO Biological Process / upregulated: production of molecular mediator of immune response (adjusted p-value=4.43e-20)
- GO Biological Process / upregulated: adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains (adjusted p-value=2.34e-19)
- GO Biological Process / upregulated: lymphocyte mediated immunity (adjusted p-value=6.67e-19)
- Reactome / upregulated: Cell Cycle Checkpoints (adjusted p-value=5.88e-17)
- GO Biological Process / upregulated: leukocyte mediated immunity (adjusted p-value=1.62e-15)
- GO Biological Process / upregulated: antimicrobial humoral response (adjusted p-value=3.87e-15)
- GO Biological Process / upregulated: antibacterial humoral response (adjusted p-value=8.11e-15)

## Interpretation caution

These pathway results should be interpreted as hypothesis-generating biological summaries of the disease-control contrast. They do not establish cell-type causality or clinical mechanism without additional validation.
