![](images/rm_logo.png)

# ReMIXTURE

Implementation of the ReMIXTURE method to quantify how (genetic or other) diversity is distributed and shared between groups, and to create pretty and intuitive plots showing this on a globe.

![](images/rmDemoVitis.png)

- Outer Circle=Total diversity
- Inner Circle=Region-unique diversity
- Line width=Overlapped diversity

ReMIXTURE is an R package, and to run it requires:

1) A symmetrical numeric matrix of pairwise sample-to-sample distances (any distance metric will do in principle), whose rownames and colnames give the region to which the sample is assigned.
2) A data.table or data.frame with column names `region`, `lon`, and `lat`, providing the position (numeric) on the globe given to each region (character).

# History

The ReMIXTURE concept was first attempted in Tripodi & Rabanus-Wallace, et al. (2021) _Global range expansion history of pepper (*Capsicum spp.*) revealed by over 10,000 genebank accessions_. PNAS. Newer versions have very significant improvements. The algorithm currently in use is not published.

# Using ReMIXTURE.

Install using `devtools::install_github("https://github.com/mtrw/ReMIXTURE")`, and follow the tutorial in ?ReMIXTURE (in the examples section).

# Future

Any questions, suggestions, feedback please email me! tim.rabanuswallace@unimelb.edu.au.
