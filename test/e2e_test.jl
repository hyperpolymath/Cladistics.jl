# SPDX-License-Identifier: MPL-2.0
# (PMPL-1.0-or-later preferred; MPL-2.0 required for Julia ecosystem)
# E2E pipeline tests for Cladistics.jl
# Tests the full phylogenetic analysis workflow: sequences → distances → trees →
# bootstrap → Newick export/import → clade identification.

using Test
using Cladistics

@testset "E2E Pipeline Tests" begin

    @testset "Full pipeline: 6-taxon phylogenetic analysis" begin
        # Representative primate-like DNA sequences (50 bp)
        sequences = [
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",  # Taxon A (reference)
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCC",  # B: 1 change
            "ATCGATCGATCAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",  # C: 1 change
            "TTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",  # D: 1 change
            "TTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCC",  # E: 2 changes
            "TTCAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",  # F: 2 changes
        ]
        taxa = ["A", "B", "C", "D", "E", "F"]

        # 1. Distance matrix — all methods
        dmat_h   = distance_matrix(sequences; method=:hamming)
        dmat_jc  = distance_matrix(sequences; method=:jc69)
        dmat_k2p = distance_matrix(sequences; method=:k2p)

        @test size(dmat_h) == (6, 6)
        @test all(dmat_h[i,i] == 0 for i in 1:6)
        @test issymmetric(dmat_h)
        @test issymmetric(dmat_jc)

        # Identical taxa have distance 0
        @test dmat_h[1,1] ≈ 0.0

        # More-different taxa have larger distance
        @test dmat_h[1,5] >= dmat_h[1,2]

        # 2. Tree inference
        tree_upgma = upgma(dmat_jc; taxa_names=taxa)
        tree_nj    = neighbor_joining(dmat_k2p; taxa_names=taxa)
        tree_mp    = maximum_parsimony(sequences; taxa_names=taxa)

        @test Set(tree_upgma.taxa) == Set(taxa)
        @test Set(tree_nj.taxa)    == Set(taxa)
        @test Set(tree_mp.taxa)    == Set(taxa)
        @test tree_upgma.method == :upgma
        @test tree_nj.method    == :nj
        @test tree_mp.method    == :parsimony

        # 3. Parsimony scoring
        cm = character_state_matrix(sequences)
        @test size(cm) == (6, length(sequences[1]))
        score = calculate_parsimony_score(tree_mp, cm)
        @test score >= 0
        @test isa(score, Int)

        # 4. Bootstrap support
        boot = bootstrap_support(sequences; replicates=20, method=:nj)
        @test all(0.0 <= v <= 1.0 for v in values(boot))

        # 5. Clade identification
        clades = identify_clades(tree_upgma, 0.0)  # threshold 0 returns all clades
        @test length(clades) >= 0

        # 6. Tree comparison
        rf = tree_distance(tree_upgma, tree_nj)
        @test rf >= 0
        @test tree_distance(tree_upgma, tree_upgma) == 0

        # 7. Root tree
        rooted = root_tree(tree_nj, "F")
        @test Set(rooted.taxa) == Set(taxa)

        # 8. Newick round-trip
        newick = tree_to_newick(rooted)
        @test endswith(newick, ";")
        reparsed = parse_newick(newick)
        @test Set(reparsed.taxa) == Set(taxa)
    end

    @testset "Error handling: invalid inputs" begin
        # Non-existent outgroup
        dmat = [0.0 0.2; 0.2 0.0]
        tree = upgma(dmat; taxa_names=["X","Y"])
        @test_throws ErrorException root_tree(tree, "NonExistent")

        # Saturated JC69 sequences produce infinite distance
        seqs_sat = ["AAAA", "TTTT"]
        dmat_sat = distance_matrix(seqs_sat; method=:jc69)
        @test isinf(dmat_sat[1,2])
    end

    @testset "Round-trip consistency: UPGMA Newick export/import" begin
        seqs = ["ATCGATCG", "ATCGATCC", "TTCGATCG", "TTCGATCC"]
        dmat = distance_matrix(seqs; method=:hamming)
        original = upgma(dmat; taxa_names=["W","X","Y","Z"])

        newick = tree_to_newick(original)
        reparsed = parse_newick(newick)

        @test Set(reparsed.taxa) == Set(original.taxa)
        @test length(reparsed.taxa) == 4
    end

end
