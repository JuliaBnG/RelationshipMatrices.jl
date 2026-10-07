"""
    _tile_size(n, colbytes) -> Int

Tile width for [`_foreach_upper_pair`](@ref). A pair of tiles holds about
4 MiB of column data (one core's share of a typical L3 slice), within
`4:128`, and is halved until there are at least `4 * nthreads()` tiles.
"""
function _tile_size(n::Integer, colbytes::Integer)
    tb = clamp(2^22 ÷ max(2colbytes, 1), 4, 128)
    ntile(tb) = (m = cld(n, tb); m * (m + 1) ÷ 2)
    while tb > 1 && ntile(tb) < 4Threads.nthreads()
        tb ÷= 2
    end
    return max(tb, 1)
end

"""
    _foreach_upper_pair(f, n, colbytes)

Call `f(i, j)` once for every pair `1 ≤ i ≤ j ≤ n`, threaded over square tiles
of the upper triangle.

`colbytes` is the size of the data one column index reads, used to pick a tile
width whose working set stays in cache. Off-diagonal tiles carry equal work
and diagonal tiles half, so even contiguous thread chunks are balanced, unlike
`@threads for j in 1:n, i in 1:j`, which leaves the last thread ≈ 2/p of the
work and streams every column from memory once per `j`.
"""
function _foreach_upper_pair(f::F, n::Integer, colbytes::Integer) where {F}
    n > 0 || return nothing
    tb = _tile_size(n, colbytes)
    m = cld(n, tb)
    tiles = [(bi, bj) for bj in 1:m for bi in 1:bj]
    Threads.@threads for t in eachindex(tiles)
        bi, bj = tiles[t]
        for j in ((bj-1)*tb+1):min(bj * tb, n), i in ((bi-1)*tb+1):min(bi * tb, j)
            f(i, j)
        end
    end
    return nothing
end
