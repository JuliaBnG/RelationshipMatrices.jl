# Segment-level IBD (roadmap Track B1–B2).
#
# Founder haplotypes carry their own labels (haplotype h of founder i is
# label 2i − 1 or 2i). Every other haplotype is a gamete of one parent: it
# starts each chromosome on one of the parent's two haplotypes and
# switches at every crossover. With the pedigree, the crossovers determine
# the founder origin of every base, as a list of segments per haplotype.
#
# Coordinates are genome-wide base positions (UInt32). Chromosome c covers
# the closed range chromosomes[c] = first:last. A crossover at base x
# keeps bases ≤ x on the current haplotype and takes bases ≥ x + 1 from
# the other one (the convention of BnGsim's `make_progeny!`). Segments are
# closed intervals [start, stop] that never cross a chromosome end;
# adjacent segments of one haplotype with the same label are merged, so
# every segment is a maximal run of one founder haplotype.

"""
    Meioses(parent, haplotype, phase, ptr, pos)
    Meioses(parent, haplotype, phase, crossovers)

The meioses that formed the non-founder haplotypes of a pedigree, one
gamete per entry `g`:

- `parent[g]`: the parent (pedigree row) that made the gamete;
- `haplotype[g]`: the haplotype it became, `2i − 1` (from the sire of
  `i`) or `2i` (from the dam);
- `phase[c, g]`: `true` if chromosome `c` starts on the parent's second
  haplotype (`2parent`), `false` for the first (`2parent − 1`);
- crossovers: genome-wide base positions, `pos[ptr[g]:ptr[g+1]-1]`,
  increasing within each chromosome; or a vector of vectors.

A crossover at base `x` keeps bases `≤ x` on the current parental
haplotype and takes bases `≥ x + 1` from the other one. BnGsim's
`MeiosisRecord` holds these fields per cohort.
"""
struct Meioses
    parent::Vector{Int32}
    haplotype::Vector{Int32}
    phase::BitMatrix
    ptr::Vector{Int}
    pos::Vector{UInt32}
    function Meioses(parent::AbstractVector{<:Integer}, haplotype::AbstractVector{<:Integer},
        phase::AbstractMatrix{Bool}, ptr::AbstractVector{<:Integer}, pos::AbstractVector{<:Integer})
        n = length(parent)
        length(haplotype) == n && size(phase, 2) == n && length(ptr) == n + 1 ||
            throw(DimensionMismatch("parent, haplotype, phase columns and ptr (n + 1) must agree"))
        ptr[1] == 1 && issorted(ptr) && ptr[end] == length(pos) + 1 ||
            throw(ArgumentError("ptr must start at 1, be non-decreasing and end at length(pos) + 1"))
        allunique(haplotype) || throw(ArgumentError("a haplotype is formed by more than one gamete"))
        new(Int32.(parent), Int32.(haplotype), BitMatrix(phase), Vector{Int}(ptr), Vector{UInt32}(pos))
    end
end

Meioses(parent, haplotype, phase, crossovers::AbstractVector{<:AbstractVector}) =
    Meioses(parent, haplotype, phase, cumsum([1; length.(crossovers)]),
        reduce(vcat, crossovers; init = UInt32[]))

Base.length(m::Meioses) = length(m.parent)

"""
    HaplotypeSegments

Founder segments of every haplotype of a pedigree, from
[`founder_segments`](@ref). Haplotype `h` (`2i − 1`, `2i` for individual
`i`) has the segments `k ∈ ptr[h]:ptr[h+1]-1`, the closed base ranges
`start[k]:stop[k]` with founder haplotype `label[k]`, in genome order,
tiling every chromosome in `chromosomes`.
"""
struct HaplotypeSegments
    chromosomes::Vector{UnitRange{UInt32}}
    ptr::Vector{Int}
    start::Vector{UInt32}
    stop::Vector{UInt32}
    label::Vector{UInt32}
end

nhaplotypes(s::HaplotypeSegments) = length(s.ptr) - 1
genome_length(s::HaplotypeSegments) = sum(length, s.chromosomes)

function Base.show(io::IO, s::HaplotypeSegments)
    print(io, "HaplotypeSegments: ", nhaplotypes(s) ÷ 2, " individuals, ", length(s.start),
        " segments on ", length(s.chromosomes), " chromosomes")
end

"""
    segments(s::HaplotypeSegments, h) -> (start, stop, label) views

The segments of haplotype `h`.
"""
function segments(s::HaplotypeSegments, h::Integer)
    r = s.ptr[h]:(s.ptr[h+1]-1)
    return view(s.start, r), view(s.stop, r), view(s.label, r)
end

"""
    founder_segments(sire, dam, meioses::Meioses, chromosomes) -> HaplotypeSegments

Founder origin of every base of every haplotype. `sire`, `dam` give the
pedigree (`0` unknown, parents before offspring); an individual is a
founder if both parents are unknown, and then its haplotypes `2i − 1`
and `2i` are founder haplotypes with those labels. Every other individual
needs both parents known and both of its haplotypes in `meioses`, each
made by the matching parent. `chromosomes` are the closed genome-wide
base ranges of the chromosomes, in order and disjoint.

Individuals are processed in pedigree order, splicing the segments of
each parent's haplotypes at the gamete's crossovers.
"""
function founder_segments(sire::AbstractVector{<:Integer}, dam::AbstractVector{<:Integer}, m::Meioses,
    chromosomes::AbstractVector{<:AbstractUnitRange})
    validate_pedigree(sire, dam)
    n = length(sire)
    chrs = [UInt32(first(r)):UInt32(last(r)) for r in chromosomes]
    all(c -> length(c) > 0, chrs) && all(k -> first(chrs[k+1]) > last(chrs[k]), 1:length(chrs)-1) ||
        throw(ArgumentError("chromosomes must be non-empty, in order and disjoint"))
    size(m.phase, 1) == length(chrs) ||
        throw(DimensionMismatch("phase has $(size(m.phase, 1)) chromosomes, chromosomes $(length(chrs))"))
    gamete = zeros(Int, 2n)                       # haplotype → gamete
    for g in eachindex(m.haplotype)
        h = m.haplotype[g]
        1 <= h <= 2n || throw(ArgumentError("gamete $g forms haplotype $h, outside 1:$(2n)"))
        gamete[h] = g
    end
    st = [UInt32[] for _ = 1:2n]
    sp = [UInt32[] for _ = 1:2n]
    lb = [UInt32[] for _ = 1:2n]
    for i = 1:n
        s, d = sire[i], dam[i]
        if s == 0 && d == 0
            for h in (2i - 1, 2i), c in chrs
                gamete[h] == 0 || throw(ArgumentError("founder $i has a gamete in meioses"))
                push!(st[h], first(c)); push!(sp[h], last(c)); push!(lb[h], UInt32(h))
            end
            continue
        end
        s > 0 && d > 0 || throw(ArgumentError("individual $i has one known parent; segments need both or none"))
        for (h, p) in ((2i - 1, s), (2i, d))
            g = gamete[h]
            g > 0 || throw(ArgumentError("haplotype $h of individual $i has no gamete in meioses"))
            m.parent[g] == p || throw(ArgumentError("gamete $g of haplotype $h comes from $(m.parent[g]), not parent $p"))
            _splice!(st[h], sp[h], lb[h], st, sp, lb, p, m, g, chrs)
        end
    end
    ptr = cumsum([1; length.(st)])
    HaplotypeSegments(chrs, ptr, reduce(vcat, st), reduce(vcat, sp), reduce(vcat, lb))
end

founder_segments(ped::DataFrame, m::Meioses, chromosomes) = founder_segments(ped.sire, ped.dam, m, chromosomes)

# Append to (ost, osp, olb) the segments of a parental haplotype
# restricted to bases a:b, merging with the last segment when contiguous,
# of the same label and on the same chromosome (starting at `cfirst`).
function _copy_range!(ost, osp, olb, pst, psp, plb, a::UInt32, b::UInt32, cfirst::UInt32)
    k = max(searchsortedlast(pst, a), 1)          # the segment containing a
    while k <= length(pst) && pst[k] <= b
        x, y = max(pst[k], a), min(psp[k], b)
        if x <= y
            if !isempty(ost) && olb[end] == plb[k] && osp[end] + 1 == x && ost[end] >= cfirst
                osp[end] = y
            else
                push!(ost, x); push!(osp, y); push!(olb, plb[k])
            end
        end
        k += 1
    end
end

function _splice!(ost, osp, olb, st, sp, lb, p, m::Meioses, g, chrs)
    xo = view(m.pos, m.ptr[g]:(m.ptr[g+1]-1))
    k = 1
    for (c, r) in enumerate(chrs)
        h = m.phase[c, g] ? 2p : 2p - 1
        a = Int(first(r))                       # Int: a crossover at the last base must not wrap
        while k <= length(xo) && xo[k] <= last(r)
            x = Int(xo[k])
            x >= first(r) || throw(ArgumentError("crossover $x of gamete $g lies outside the chromosomes or is out of order"))
            if x >= a
                _copy_range!(ost, osp, olb, st[h], sp[h], lb[h], UInt32(a), UInt32(x), first(r))
                a = x + 1
            end
            h = h == 2p ? 2p - 1 : 2p
            k += 1
        end
        a <= last(r) && _copy_range!(ost, osp, olb, st[h], sp[h], lb[h], UInt32(a), last(r), first(r))
    end
    k == length(xo) + 1 || throw(ArgumentError("crossovers of gamete $g lie outside the chromosomes or are out of order"))
end

# Bases shared IBD by haplotypes a and b, counting only shared runs of
# at least `min_len` bases. A shared run is a maximal stretch over which
# the two haplotypes carry the same founder haplotype, merged across
# founder labels: when both switch founder at the same base, they
# inherited that junction from a common ancestor, so the stretch is one
# segment identical by descent from that ancestor. Runs do not cross a
# chromosome end. Segments of each haplotype are in genome order and
# never cross a chromosome end, so one two-pointer pass covers the genome.
function _shared(s::HaplotypeSegments, a::Integer, b::Integer, min_len::Integer)
    i, ie = s.ptr[a], s.ptr[a+1] - 1
    j, je = s.ptr[b], s.ptr[b+1] - 1
    tot = 0
    run = 0                                       # length of the current shared run
    last_y = UInt32(0)                            # its last base
    c = 1                                         # chromosome of the current position
    chrs = s.chromosomes
    @inbounds while i <= ie && j <= je
        x = max(s.start[i], s.start[j])
        y = min(s.stop[i], s.stop[j])
        if x <= y && s.label[i] == s.label[j]
            while x > last(chrs[c])
                c += 1
            end
            if run > 0 && x == last_y + 1 && x != first(chrs[c])
                run += Int(y - x) + 1             # contiguous: the same run
            else
                run >= min_len && (tot += run)
                run = Int(y - x) + 1
            end
            last_y = y
        end
        s.stop[i] < s.stop[j] ? (i += 1) : s.stop[i] > s.stop[j] ? (j += 1) : (i += 1; j += 1)
    end
    run >= min_len && (tot += run)
    return tot
end

"""
    irm(segs::HaplotypeSegments; ids = all, min_len = 0, T = Float64) -> Matrix{T}

Realized IBD relationship matrix of the individuals `ids` from founder
segments, on the scale of the numerator relationship matrix:

```math
R_{ij} = \\frac{1}{2G}\\sum_{a \\in \\{i_1, i_2\\}}\\sum_{b \\in \\{j_1, j_2\\}} S_{ab},
```

``S_{ab}`` the number of bases at which haplotypes ``a`` and ``b`` carry
the same founder haplotype and ``G`` the genome length in bases. The
diagonal is ``1 + F_i``. This is the locus-level [`irm`](@ref) with a
locus at every base, without sampling loci. With `min_len`, only shared
runs of at least `min_len` bases count (recent IBD). A shared run is a
maximal stretch identical by descent, merged across founder labels (two
haplotypes that switch founder at the same base inherited that junction
from a common ancestor) and ending at chromosome ends.
"""
function irm(segs::HaplotypeSegments; ids::AbstractVector{<:Integer} = 1:(nhaplotypes(segs)÷2),
    min_len::Integer = 0, T::Type{<:AbstractFloat} = Float64)
    nid = nhaplotypes(segs) ÷ 2
    all(i -> 1 <= i <= nid, ids) || throw(ArgumentError("ids must lie in 1:$nid"))
    n = length(ids)
    G = genome_length(segs)
    R = Matrix{T}(undef, n, n)
    perhap = 12 * max(1, length(segs.start) ÷ max(1, nhaplotypes(segs)))
    _foreach_upper_pair(n, 2perhap) do x, y
        i, j = ids[x], ids[y]
        tot = 0
        for a in (2i - 1, 2i), b in (2j - 1, 2j)
            tot += a == b ? G : _shared(segs, a, b, min_len)
        end
        R[x, y] = R[y, x] = T(tot / 2G)
    end
    return R
end

"""
    segment_inbreeding(segs::HaplotypeSegments; ids = all, min_len = 0) -> Vector{Float64}

Realized inbreeding coefficient of each individual: the proportion of the
genome at which its two haplotypes carry the same founder haplotype
(shared runs of at least `min_len` bases, as in [`irm`](@ref) on
segments: the autozygous segments that runs of homozygosity estimate). Equals the diagonal of
[`irm`](@ref) on segments minus one.
"""
function segment_inbreeding(segs::HaplotypeSegments; ids::AbstractVector{<:Integer} = 1:(nhaplotypes(segs)÷2),
    min_len::Integer = 0)
    G = genome_length(segs)
    [_shared(segs, 2i - 1, 2i, min_len) / G for i in ids]
end
