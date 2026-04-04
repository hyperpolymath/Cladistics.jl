# SPDX-License-Identifier: MPL-2.0
# (PMPL-1.0-or-later preferred; MPL-2.0 required for Julia ecosystem)
# BenchmarkTools benchmarks for Cladistics.jl
# Measures distance matrix construction, tree inference, and bootstrap at
# small/medium/large scales.

using BenchmarkTools
using Cladistics

# ── Helpers ───────────────────────────────────────────────────────────────────

rand_seq(len::Int) = join(rand(['A','T','C','G'], len))

function mutate_seq(seq::String, k::Int)
    s = collect(seq)
    for p in randperm(length(s))[1:min(k, length(s))]
        s[p] = rand(['A','T','C','G'])
    end
    join(s)
end

function make_seqs(n_taxa::Int, seq_len::Int, max_mut::Int)
    base = rand_seq(seq_len)
    [mutate_seq(base, rand(0:max_mut)) for _ in 1:n_taxa]
end

# Small: 5 taxa, 50 bp
seqs_small  = make_seqs(5,  50,  3)
# Medium: 20 taxa, 200 bp
seqs_medium = make_seqs(20, 200, 8)
# Large: 50 taxa, 500 bp
seqs_large  = make_seqs(50, 500, 20)

dmat_small  = distance_matrix(seqs_small;  method=:hamming)
dmat_medium = distance_matrix(seqs_medium; method=:jc69)
dmat_large  = distance_matrix(seqs_large;  method=:k2p)

# ── Distance matrix benchmarks ────────────────────────────────────────────────

println("=== distance_matrix :hamming (small: 5×50) ===")
@benchmark distance_matrix($seqs_small; method=:hamming)

println("=== distance_matrix :jc69 (medium: 20×200) ===")
@benchmark distance_matrix($seqs_medium; method=:jc69)

println("=== distance_matrix :k2p (large: 50×500) ===")
@benchmark distance_matrix($seqs_large; method=:k2p)

# ── UPGMA tree benchmarks ─────────────────────────────────────────────────────

println("=== upgma (small: 5 taxa) ===")
@benchmark upgma($dmat_small)

println("=== upgma (medium: 20 taxa) ===")
@benchmark upgma($dmat_medium)

println("=== upgma (large: 50 taxa) ===")
@benchmark upgma($dmat_large)

# ── NJ tree benchmarks ────────────────────────────────────────────────────────

println("=== neighbor_joining (small: 5 taxa) ===")
@benchmark neighbor_joining($dmat_small)

println("=== neighbor_joining (medium: 20 taxa) ===")
@benchmark neighbor_joining($dmat_medium)

println("=== neighbor_joining (large: 50 taxa) ===")
@benchmark neighbor_joining($dmat_large)
