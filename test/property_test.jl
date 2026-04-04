# SPDX-License-Identifier: MPL-2.0
# (PMPL-1.0-or-later preferred; MPL-2.0 required for Julia ecosystem)
# Property-based tests for Cladistics.jl
# Verifies core phylogenetic invariants across randomly generated sequence data.

using Test
using Cladistics

@testset "Property-Based Tests" begin

    # Helper: generate a random DNA sequence of given length
    function rand_seq(len::Int)
        join(rand(['A','T','C','G'], len))
    end

    # Helper: introduce k random substitutions into a sequence
    function mutate_seq(seq::String, k::Int)
        s = collect(seq)
        positions = randperm(length(s))[1:min(k, length(s))]
        for p in positions
            s[p] = rand(['A','T','C','G'])
        end
        join(s)
    end

    @testset "Invariant: distance matrix is symmetric and zero-diagonal" begin
        for _ in 1:50
            n_taxa = rand(3:6)
            seq_len = rand(8:20)
            base = rand_seq(seq_len)
            seqs = [mutate_seq(base, rand(0:3)) for _ in 1:n_taxa]
            for method in [:hamming, :p_distance]
                dmat = distance_matrix(seqs; method=method)
                @test size(dmat) == (n_taxa, n_taxa)
                @test dmat ≈ dmat'  atol=1e-12
                @test all(dmat[i,i] ≈ 0.0 for i in 1:n_taxa)
            end
        end
    end

    @testset "Invariant: distance matrix entries are non-negative" begin
        for _ in 1:50
            n_taxa = rand(2:5)
            seqs = [rand_seq(rand(6:15)) for _ in 1:n_taxa]
            # Ensure uniform length
            min_len = minimum(length.(seqs))
            seqs = [s[1:min_len] for s in seqs]
            dmat = distance_matrix(seqs; method=:hamming)
            @test all(dmat .>= 0.0)
        end
    end

    @testset "Invariant: UPGMA tree preserves all taxa" begin
        for _ in 1:50
            n_taxa = rand(3:7)
            seq_len = rand(8:16)
            base = rand_seq(seq_len)
            seqs = [mutate_seq(base, rand(0:3)) for _ in 1:n_taxa]
            dmat = distance_matrix(seqs; method=:hamming)
            taxa_names = ["T$i" for i in 1:n_taxa]
            tree = upgma(dmat; taxa_names=taxa_names)
            @test Set(tree.taxa) == Set(taxa_names)
            @test length(tree.taxa) == n_taxa
        end
    end

    @testset "Invariant: NJ tree preserves all taxa" begin
        for _ in 1:50
            n_taxa = rand(3:7)
            seq_len = rand(8:16)
            base = rand_seq(seq_len)
            seqs = [mutate_seq(base, rand(0:3)) for _ in 1:n_taxa]
            dmat = distance_matrix(seqs; method=:p_distance)
            taxa_names = ["S$i" for i in 1:n_taxa]
            tree = neighbor_joining(dmat; taxa_names=taxa_names)
            @test Set(tree.taxa) == Set(taxa_names)
        end
    end

    @testset "Invariant: parsimony score is non-negative" begin
        for _ in 1:50
            n_taxa = rand(3:6)
            seq_len = rand(6:12)
            base = rand_seq(seq_len)
            seqs = [mutate_seq(base, rand(0:3)) for _ in 1:n_taxa]
            dmat = distance_matrix(seqs; method=:hamming)
            taxa_names = ["R$i" for i in 1:n_taxa]
            tree = upgma(dmat; taxa_names=taxa_names)
            cm = character_state_matrix(seqs)
            score = calculate_parsimony_score(tree, cm)
            @test score >= 0
        end
    end

    @testset "Invariant: bootstrap support values in [0, 1]" begin
        for _ in 1:50
            n_taxa = rand(3:5)
            seq_len = rand(10:20)
            base = rand_seq(seq_len)
            seqs = [mutate_seq(base, rand(0:3)) for _ in 1:n_taxa]
            boot = bootstrap_support(seqs; replicates=10, method=:upgma)
            for (_, v) in boot
                @test 0.0 <= v <= 1.0
            end
        end
    end

end
