# Bachelor's Thesis: Improving the Tumor Specificity of CAR-T Cell Therapy for Neuroblastoma by Finding Suitable Combinations of Target Antigens

This repository contains the code related to my bachelor's thesis where I applied a genetic algorithm for combinatorial
optimization and random forest feature selection for binary classification in a bioinformatics context.

## Abstract

In recent years, CAR-T cell therapy has led to remarkable clinical success in treating hematological malignancies like
leukemia. In order to achieve similar results in solid tumors like neuroblastoma, however, additional barriers must be
overcome. One of these barriers is the lack of truly tumor-specific antigens, which can lead to severe on-target,
off-tumor toxicities. Therefore, novel methods for increasing target specificity are in development. This includes
AND-gate CARs, which target multiple antigens and only activate if all of them are expressed simultaneously on a cell.

In this thesis, we use RNA gene expression data to determine sets of genes that could be promising candidates for
combinatorial targeting based on their disparity in RNA expression levels between cancerous and healthy tissues. We
frame the problem as (a) a combinatorial optimization problem, and (b) a binary classification problem. We formally
define an objective function and evaluate the suitability of (a) the genetic algorithm (GA) metaheuristic, and (b)
random forest (RF) feature selection for finding high-quality solutions. Additionally, we verify the plausibility of our
results using reference antigens studied in previous work.

We find that although both approaches yield good results compared to a random baseline, the GA solutions consistently
outperform those selected by the RF. Furthermore, we observe that a substantial separation of neuroblastoma samples from
healthy tissue samples can be achieved using sets of only three genes, although solution quality continues to improve
with an increasing number of selected genes. Finally, we find that the frequent occurrence of a majority of reference
antigens in our results corroborates the initial hypothesis that RNA expression data combined with the employed
objective function could be a useful tool for finding combinations of target antigens for use _in vivo_.