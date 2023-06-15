# Multi-Ancestry Meta-Analysis and Fine-Mapping in Parkinson's Disease 

`GP2 ❤️ Open Science 😍`

 - **Project:** Multi-ancestry Meta-Analysis and Fine-Mapping in PD
 - **Date Last Updated:** June 2023 
    - **Update Description:** Updated README, relocated repo

---

## Summary
This is the online repository for the manuscript titled ***"Multi-ancestry genome-wide meta-analysis in Parkinson’s disease"***. This study is the first large-scale multi-ancestry meta-analysis of PD GWASs, incorporating data from 49,049 PD cases, 18,785 proxy cases, and 2,458,063 controls including individuals of European, East Asian, Latin American, and African ancestry. 

### Data Statement 
Data used include:
- European: **Identification of novel risk loci, causal insights, and heritable risk for Parkinson's disease: a meta-analysis of genome-wide association studies.** *Nalls et al. 2019* ([PubMed](https://pubmed.ncbi.nlm.nih.gov/31701892/))
- Latin American: **LARGE‐PD: Examining the genetics of Parkinson's disease in Latin America.** *Loesch et al. 2021* ([PubMed](https://pubmed.ncbi.nlm.nih.gov/34227697/))
- Asian: **Identification of Risk Loci for Parkinson Disease in Asians and Comparison of Risk Between Asians and Europeans: A Genome-Wide Association Study.** *Foo et al. 2020* ([PubMed](https://pubmed.ncbi.nlm.nih.gov/32310270/))
- 23andMe GWAS summary statistics (available via collaboration with 23andMe).

### Helpful Links 
- [medRxiv pre-print (August 2022)](https://doi.org/10.1101/2022.08.04.22278432)
- [GP2 Website](https://gp2.org/)
    - [GP2 Cohort Dashboard](https://gp2.org/cohort-dashboard-advanced/)
- [Introduction to GP2](https://movementdisorders.onlinelibrary.wiley.com/doi/10.1002/mds.28494)
    - [Other GP2 Manuscripts (PubMed)](https://pubmed.ncbi.nlm.nih.gov/?term=%22global+parkinson%27s+genetics+program%22)


# Workflow Diagram 
![Figure 1](https://github.com/GP2code/GP2-Multiancestry-metaGWAS/assets/44064705/13cc9d34-3658-4f75-8387-15b9a8235ed3)

---
# Analysis Notebooks
* Languages: Python and bash

| **Notebooks** |                                                    **Description**                                                   |
|:----------------:|:--------------------------------------------------------------------------------------------------------------------:|
|        PD_MAMA.PT1.Meta-Analysis.md    | Meta-analysis |
|        PD_MAMA.PT2.Finemap.md  | Fine-mapping |

---
### Figures and Supplemental Figures

![Figure 2](https://github.com/GP2code/GP2-Multiancestry-metaGWAS/assets/44064705/b98bdef9-c4f8-4668-bfb6-5ea21eeac74e)

Manhattan plots and upset plot of PD MAMA results. A: random-effect; B: MR-MEGA. Orange dotted line indicates the Bonferroni adjusted significant threshold of P < 5 x 10<sup>-9</sup>. Gray dotted line indicates the truncation line, where all -log10P values greater than 40 were truncated to 40 for visual clarity. Novel loci are highlighted in red and annotated with the nearest protein coding gene. C: Heterogeneity upset plot of the top hits in novel loci. The top bar plot illustrates heterogeneity with dark blue indicating ancestry heterogeneity proportion and light blue indicating other residual heterogeneity proportion. The bottom plot shows the subcohort level beta values with blue indicating positive and redindicating negative effect directions. 3 variants with greater than 30% I2 total heterogeneity were only identified in the MR-MEGA meta-analysis method, while little to no heterogeneity is observed in loci identified in random-effect. D: Heterogeneity upset plot of the top variant per MR-MEGA identified locus that had moderate to high heterogeneity (I<sup>2</sup> > 30). Variants in novel loci are annotated with *


---

# Software 

|               Software              |      Version(s)     |                       Resource URL                       |       RRID      |                                               Notes                                               |
|:-----------------------------------:|:-------------------:|:--------------------------------------------------------:|:---------------:|:-------------------------------------------------------------------------------------------------:|
|                PLINK                | 1.9 |            http://www.nitrc.org/projects/plink           | RRID:SCR_001757 |                                     used for fixed/random effect meta-analyses                                     |
|                MR-MEGA                | v.0.2 |            https://www.geenivaramu.ee/en/tools/mr-mega           | N/A |                                     used for meta-analysis and fine-mapping                                     |
|     Python Programming Language     |     3.8 and 3.9     |                  http://www.python.org/                  | RRID:SCR_008394 | pandas; numpy; seaborn; matplotlib; statsmodel; used for general data wrangling/plotting/analyses |
| R Project for Statistical Computing |         4.2         |                 http://www.r-project.org/                | RRID:SCR_001905 |   tidyverse; dplyr; tidyr; ggplot; data.table; used for general data wrangling/plotting/analyses  |
