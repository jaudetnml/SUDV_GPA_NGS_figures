# Filtering and plotting of NGS results for the adaptation of SUDV to guinea pigs

This repository accompanies the publication:

**Submitted, no actual ref yet**

## Analysis details

The data were analyzed as described in the publication using Geneious.

This repository uses the exported depths and vcf files to generate a SNP heat map
and line graphs used in the publication.

In order to be accepted considered for plotting, a mutation must be in a region
with at least 1000 reads and reach a frequency of at least 5%.
Lower frequencies are plotted as 0.
Regions with lower read counts are considered missing data.

Important: In the SNP_effects.csv file, codons that change to stop codons have a Protein effect of X to NA mutation.
This is due to Geneious not exporting the codon change when a stop codon is introduced.

## Running the analysis

To generate the figures from the data, switch to the R folder and run the script:

```
cd R
Rscript --vanilla --no-save ./PlotSNP_Pub.r
``` 

## Dependencies

The following packages were used to run this script:

 * tidyverse  v2.0.0
 * cowplot  v1.1.3
 * patchwork  v1.2.0
 * bioseq (from bioconductor) v0.1.4
 * ggrepel  v0.9.5
 * writexl  v1.5.0

--------

# Legal

Copyright Government of Canada 2024

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed uder the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language gorverning permissions and limitations uder the License.

--------

# CHANGELOG

20241011 : Initial commit of code and data for publication. Commit tagged as v1.0.0