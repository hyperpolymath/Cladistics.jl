[![Sponsor](https://img.shields.io/badge/Sponsor-%E2%9D%A4-pink?logo=github)](https://github.com/sponsors/hyperpolymath)

// SPDX-License-Identifier: MPL-2.0
// SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>

= Cladistics.jl
:toc: preamble
:icons: font

image:https://img.shields.io/badge/OpenSSF-Best_Practices-green?logo=opensourcesecurity[OpenSSF Best Practices,link="https://www.bestpractices.dev/en/projects/new?repo_url=https://github.com/hyperpolymath/Cladistics.jl"]
image:https://img.shields.io/badge/License-PMPL--1.0-blue.svg[License: PMPL-1.0,link="https://github.com/hyperpolymath/palimpsest-license"]
image:https://api.thegreenwebfoundation.org/greencheckimage/github.com[Green Web,link="https://www.thegreenwebfoundation.org/green-web-check/?url=github.com"]
image:https://img.shields.io/badge/Project-Topology-9558B2[Project Topology,link="TOPOLOGY.md"]
image:https://img.shields.io/badge/Completion-75%25-yellow[Completion Status,link="TOPOLOGY.md"]
image:https://img.shields.io/badge/license-PMPL--1.0--or--later-blue.svg[License,link="LICENSE"]
image:https://img.shields.io/badge/julia-1.6+-purple.svg[Julia,link="https://julialang.org"]

A Julia package for phylogenetic analysis and cladistics — the study of
evolutionary relationships among organisms.

== Overview

Cladistics is a method of biological classification that groups organisms based
on their evolutionary ancestry and shared derived characteristics (synapomorphies).
This package provides computational tools for reconstructing and analyzing
phylogenetic trees from molecular sequence data and morphological characters.

=== What is Cladistics?

* *Clade*: A monophyletic group consisting of an ancestor and all its descendants
* *Synapomorphy*: A shared derived characteristic that defines a clade
* *Phylogenetic Tree*: A branching diagram showing evolutionary relationships
* *Parsimony*: The principle that the simplest evolutionary explanation (fewest changes) is preferred
* *Bootstrap Analysis*: Statistical method to assess confidence in tree topology

== Features

=== Distance-Based Methods

* *Multiple Distance Metrics*: Hamming, p-distance, Jukes-Cantor 1969, Kimura 2-parameter
* *UPGMA*: Unweighted Pair Group Method with Arithmetic Mean (assumes molecular clock)
* *Neighbor-Joining*: Does not assume molecular clock, handles rate variation

=== Character-Based Methods

* *Maximum Parsimony*: Find trees requiring fewest evolutionary changes
* *Fitch Algorithm*: Efficient parsimony score calculation
* *Parsimony-Informative Sites*: Identify characters useful for phylogenetic inference

=== Tree Analysis

* *Bootstrap Support*: Assess confidence in tree topology (1000+ replicates)
* *Clade Identification*: Extract well-supported monophyletic groups
* *Robinson-Foulds Distance*: Compare tree topologies quantitatively
* *Tree Rooting*: Root unrooted trees using outgroup taxa
* *Newick Format*: Export trees in standard phylogenetic format

== Installation

[source,julia]
----
using Pkg
Pkg.add("Cladistics")
----

== Quick Start

[source,julia]
----
using Cladistics
using Plots

# 1. DNA sequence data (aligned)
sequences = [
    "ATCGATCGATCG",  # Species A
    "ATCGATCGATCG",  # Species B (same as A)
    "ATCGTTCGATCG",  # Species C
    "TTCGATCGATCG",  # Species D
    "TTCGTTCGATCC"   # Species E (outgroup)
]
taxa_names = ["Species_A", "Species_B", "Species_C", "Species_D", "Species_E"]

# 2. Calculate evolutionary distances
dmat = distance_matrix(sequences, method=:jc69)

# 3. Build phylogenetic tree using UPGMA
upgma_tree = upgma(dmat, taxa_names=taxa_names)
newick = tree_to_newick(upgma_tree)

# 4. Build tree using Neighbor-Joining
nj_tree = neighbor_joining(dmat, taxa_names=taxa_names)

# 5. Compare tree topologies
rf_distance = tree_distance(upgma_tree, nj_tree)

# 6. Bootstrap analysis for confidence
support = bootstrap_support(sequences, replicates=1000, method=:nj)

# 7. Identify well-supported clades
clades = identify_clades(upgma_tree, 0.95)
----

== Distance Methods Comparison

[source,julia]
----
using Cladistics
sequences = ["ATCG", "ATCG", "TTCG", "TTCC"]

hamming = distance_matrix(sequences, method=:hamming)   # Simple counting
p_dist  = distance_matrix(sequences, method=:p_distance) # Proportion of differences
jc69    = distance_matrix(sequences, method=:jc69)       # Jukes-Cantor
k2p     = distance_matrix(sequences, method=:k2p)        # Kimura 2-parameter
----

== Tree Construction Methods

=== UPGMA: Simple but assumes molecular clock

[source,julia]
----
dmat = [0.0 0.2 0.4; 0.2 0.0 0.3; 0.4 0.3 0.0]
tree = upgma(dmat, taxa_names=["Human", "Chimp", "Gorilla"])
# Good for: recent divergences, molecular clock valid
# Not good for: ancient divergences, variable rates
----

=== Neighbor-Joining: No molecular clock assumption

[source,julia]
----
dmat = [0.0 0.5 0.8 1.0; 0.5 0.0 0.6 0.9; 0.8 0.6 0.0 0.4; 1.0 0.9 0.4 0.0]
tree = neighbor_joining(dmat, taxa_names=["A", "B", "C", "D"])
# Good for: variable evolutionary rates, large datasets
----

== Bootstrap Analysis

[source,julia]
----
using Cladistics, Random
Random.seed!(42)
sequences = [
    "ATCGATCGATCGATCG",
    "ATCGATCGATCGATCG",
    "ATCGTTCGATCGATCG",
    "TTCGATCGATCGATCG",
    "TTCGTTCGATCGTTCC"
]
# 1000 bootstrap replicates — resamples alignment columns with replacement
support = bootstrap_support(sequences, replicates=1000, method=:nj)
# Interpret: >95% Strong, 70-95% Moderate, <70% Weak
----

== Maximum Parsimony

[source,julia]
----
using Cladistics
alignment = ["ATCG", "ATCG", "TTCG", "TTCC"]
char_matrix = character_state_matrix(alignment)
informative_sites = parsimony_informative_sites(char_matrix)
dmat = distance_matrix(alignment, method=:hamming)
tree = upgma(dmat, taxa_names=["Tax1", "Tax2", "Tax3", "Tax4"])
score = calculate_parsimony_score(tree, char_matrix)
# Lower score = more parsimonious (preferred)
----

== Real-World Example: Primate Phylogeny

[source,julia]
----
using Cladistics
primate_sequences = [
    "ATCGATCGATCGATCGATCG",  # Human
    "ATCGATCGATCGATCGATCG",  # Chimp
    "ATCGATCGTTCGATCGATCG",  # Gorilla
    "ATCGTTCGTTCGATCGTTCG",  # Orangutan
    "TTCGTTCGTTCGTTCGTTCG"   # Lemur (outgroup)
]
taxa = ["Human", "Chimp", "Gorilla", "Orangutan", "Lemur"]
dmat = distance_matrix(primate_sequences, method=:k2p)
tree = neighbor_joining(dmat, taxa_names=taxa)
rooted_tree = root_tree(tree, "Lemur")
support = bootstrap_support(primate_sequences, replicates=100, method=:nj)
newick = tree_to_newick(rooted_tree)
----

== Key Concepts

=== Molecular Clock

* Assumption: constant evolutionary rate across lineages
* When valid: recent species, similar generation times
* Methods: UPGMA assumes molecular clock
* When violated: use Neighbor-Joining

=== Bootstrap Support Interpretation

* `>95%`: Publish with confidence
* `70–95%`: Mention uncertainty
* `<70%`: Weak support, collect more data

=== Parsimony vs Distance

* *Parsimony*: Find tree requiring fewest character changes; NP-hard but theoretically clear
* *Distance*: Fast, scales to large datasets; loses some information from pairwise comparisons

== References

=== Classic Papers

* Felsenstein, J. (1985). "Confidence limits on phylogenies: An approach using the bootstrap." _Evolution_, 39(4), 783–791.
* Saitou, N., & Nei, M. (1987). "The neighbor-joining method." _Molecular Biology and Evolution_, 4(4), 406–425.
* Fitch, W. M. (1971). "Toward defining the course of evolution." _Systematic Zoology_, 20(4), 406–416.

=== Textbooks

* Felsenstein, J. (2004). _Inferring Phylogenies_. Sinauer Associates.
* Lemey, P. et al. (2009). _The Phylogenetic Handbook_. Cambridge University Press.
* Hall, B. G. (2011). _Phylogenetic Trees Made Easy_. Sinauer Associates.

=== Evolutionary Models

* Jukes, T. H., & Cantor, C. R. (1969). "Evolution of protein molecules." _Mammalian Protein Metabolism_, pp. 21–132.
* Kimura, M. (1980). "A simple method for estimating evolutionary rates." _Journal of Molecular Evolution_, 16(2), 111–120.

== Citation

[source,bibtex]
----
@software{cladistics_jl,
  author = {Jewell, Jonathan D.A.},
  title = {Cladistics.jl: Phylogenetic Analysis in Julia},
  year = {2026},
  url = {https://github.com/hyperpolymath/Cladistics.jl}
}
----

== Related Projects

* https://github.com/BioJulia[BioJulia] — Broader bioinformatics ecosystem
* https://github.com/crsl4/PhyloNetworks.jl[PhyloNetworks.jl] — Phylogenetic networks
* https://github.com/richardreeve/Phylo.jl[Phylo.jl] — Alternative phylogenetics package

== External Visualization Tools

* http://tree.bio.ed.ac.uk/software/figtree/[FigTree]
* https://itol.embl.de/[iTOL]
* http://etetoolkit.org/[ETE Toolkit]

== Contributing

See link:CONTRIBUTING.md[CONTRIBUTING.md] for guidelines.

Wondering how this works? See link:EXPLAINME.adoc[].

== License

SPDX-License-Identifier: MPL-2.0 +
See link:LICENSE[LICENSE].
