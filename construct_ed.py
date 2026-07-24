#!/usr/bin/env python3
"""
ED string construction from synteny-block FASTA families.

Input FASTA header example:
    >5653|GCF_000005845.2_ASM584v2_genomic.fna|NC_000913.3:1-5593(-)

The number before the first "|" is used as the synteny-block identifier. All
sequences sharing an identifier form one block and are collapsed into a single
elastic-degenerate (ED) string.

The outputs of the upstream extraction step are consumed here directly.
Reference: https://github.com/kiriltemelkov/synteny-block-sequence-extraction

Output: one ED string per synteny block, in FASTA format.

ED format example:
    >A_block_id
    {A,C,}GAAT{AT,A}ATT

Meaning:
    - deterministic segments are written as plain strings,
    - nondeterministic segments are written as {variant1,variant2,...},
    - the empty word is represented by an empty variant, e.g. {A,C,}.

    
Algorithm:

For each block:

    1. Remove duplicate whole sequences.
    2. If a single distinct sequence remains, emit it
    3. Otherwise anchor and decompose recursively:
       a. Build a generalized suffix array using SA-IS over the sequences, separated by pairwise-distinct separator symbols and terminated by one unique final sentinel.
       b. Compute the LCP array using Kasai's algorithm.
       c. Find all multiMUMs: maximal exact matches that occur in every sequence and exactly once in each sequence.
       d. Select a maximum-weight collinear, non-overlapping anchor chain.
       e. Emit the anchors as deterministic segments and recursively process the gaps between them. Gaps without a common anchor collapse into one nondeterministic segment.

    
Only the Python standard library is used.
"""

from __future__ import annotations

import argparse
import collections
import concurrent.futures
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple



# Data structures

@dataclass(frozen=True)
class FastaRecord:
    header: str
    block_id: str
    seq: str


@dataclass(frozen=True)
class Anchor:
    """
    A common exact anchor.

    positions[i] is the start of the anchor in lane i, length is its length and text is the (identical) anchor sequence. Because anchors are multiMUMs there
    is exactly one occurrence per lane, so a single position tuple describes it.
    """
    positions: Tuple[int, ...]
    length: int
    text: str



# FASTA parsing and writing

def parse_block_id(header: str) -> str:
    """
    Extract the synteny-block identifier from a FASTA header.

    Example:
        5653|GCF_000005845.2_ASM584v2_genomic.fna|NC_000913.3:1-5593(-) -> 5653
    """
    return header.split("|", 1)[0].strip()


def read_fasta(path: str) -> List[FastaRecord]:
    records: List[FastaRecord] = []
    header: Optional[str] = None
    seq_parts: List[str] = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line[0] == ">":
                if header is not None:
                    seq = "".join(seq_parts).upper()
                    records.append(FastaRecord(header, parse_block_id(header), seq))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())

    if header is not None:
        seq = "".join(seq_parts).upper()
        records.append(FastaRecord(header, parse_block_id(header), seq))

    return records


def write_fasta(records: List[Tuple[str, str]], path: str) -> None:
    """
    Write ED records. Each ED string is emitted on a single line. 
    """
    with open(path, "w", encoding="utf-8") as out:
        for header, seq in records:
            out.write(">")
            out.write(header)
            out.write("\n")
            out.write(seq)
            out.write("\n")



# Generalized text encoding

def build_generalized(lanes: Sequence[str]) -> Tuple[List[int], List[int], List[int], int]:
    """
    Encode lanes as one integer text with pairwise-distinct sentinels.

    Codes:
        0                 : final sentinel (unique global minimum, required by SA-IS)
        1 .. A            : the A distinct input characters
        A+1 .. A+k        : one distinct separator per lane

    Distinct separators guarantee that no common substring crosses a lane boundary, which is exactly the generalized-suffix-array property we need.

    Returns:
        text      : integer text, terminated by the sentinel
        owner     : owner[p] = lane index, or -1 for separators / sentinel
        local_pos : local_pos[p] = position inside the lane, or -1
        alphabet_size
    """
    alphabet = sorted({ch for lane in lanes for ch in lane})
    char_code = {ch: i + 1 for i, ch in enumerate(alphabet)}
    a = len(alphabet)
    k = len(lanes)

    text: List[int] = []
    owner: List[int] = []
    local: List[int] = []

    for sid, lane in enumerate(lanes):
        for j, ch in enumerate(lane):
            text.append(char_code[ch])
            owner.append(sid)
            local.append(j)
        text.append(a + 1 + sid)  # lane separator (distinct per lane)
        owner.append(-1)
        local.append(-1)

    text.append(0)  # final sentinel
    owner.append(-1)
    local.append(-1)

    alphabet_size = a + k + 1
    return text, owner, local, alphabet_size



# Suffix array (SA-IS, linear time) and LCP (Kasai, linear time)

def sa_is(s: List[int], k: int) -> List[int]:
    """
    Linear-time suffix array via induced sorting (Nong, Zhang & Chan, 2009).

    `s` must be a list of integers in [0, k) whose last element is 0 and is the unique global minimum (the sentinel).
    """
    n = len(s)
    if n == 0:
        return []
    if n == 1:
        return [0]

    # Type map: True == S-type, False == L-type.
    t = bytearray(n)
    t[n - 1] = 1
    for i in range(n - 2, -1, -1):
        if s[i] < s[i + 1] or (s[i] == s[i + 1] and t[i + 1]):
            t[i] = 1

    is_lms = bytearray(n)
    for i in range(1, n):
        if t[i] and not t[i - 1]:
            is_lms[i] = 1
    lms_text = [i for i in range(1, n) if is_lms[i]]

    def bucket_bounds(at_end: bool) -> List[int]:
        counts = [0] * k
        for c in s:
            counts[c] += 1
        bounds = [0] * k
        acc = 0
        for i in range(k):
            acc += counts[i]
            bounds[i] = acc if at_end else acc - counts[i]
        return bounds

    def induce(sa: List[int]) -> None:
        # Induce L-type suffixes left to right from bucket starts.
        bkt = bucket_bounds(False)
        for i in range(n):
            j = sa[i] - 1
            if j >= 0 and not t[j]:
                c = s[j]
                sa[bkt[c]] = j
                bkt[c] += 1
        # Induce S-type suffixes right to left from bucket ends.
        bkt = bucket_bounds(True)
        for i in range(n - 1, -1, -1):
            j = sa[i] - 1
            if j >= 0 and t[j]:
                c = s[j]
                bkt[c] -= 1
                sa[bkt[c]] = j

    # Step 1: place LMS suffixes (any order within a bucket) and induce.
    sa = [-1] * n
    bkt = bucket_bounds(True)
    for p in reversed(lms_text):
        c = s[p]
        bkt[c] -= 1
        sa[bkt[c]] = p
    induce(sa)

    # Step 2: name LMS substrings in sorted order.
    lms_sorted = [p for p in sa if is_lms[p]]
    n1 = len(lms_sorted)
    name_of = [-1] * n
    name = 0
    prev = -1
    for p in lms_sorted:
        if prev == -1:
            name_of[p] = 0
            prev = p
            continue
        diff = False
        d = 0
        while True:
            pp = p + d
            qq = prev + d
            if pp >= n or qq >= n:
                diff = True
                break
            if s[pp] != s[qq] or t[pp] != t[qq]:
                diff = True
                break
            a_end = d > 0 and is_lms[pp]
            b_end = d > 0 and is_lms[qq]
            if a_end and b_end:
                break
            if a_end != b_end:
                diff = True
                break
            d += 1
        if diff:
            name += 1
        name_of[p] = name
        prev = p
    num_names = name + 1

    # Step 3: recurse on the reduced string if names are not already unique.
    reduced = [name_of[p] for p in lms_text]
    if num_names == n1:
        sa1 = [0] * n1
        for i, nm in enumerate(reduced):
            sa1[nm] = i
    else:
        sa1 = sa_is(reduced, num_names)

    # Step 4: place LMS suffixes in their true sorted order and induce.
    sa = [-1] * n
    bkt = bucket_bounds(True)
    for i in range(n1 - 1, -1, -1):
        p = lms_text[sa1[i]]
        c = s[p]
        bkt[c] -= 1
        sa[bkt[c]] = p
    induce(sa)
    return sa


def kasai_lcp(s: List[int], sa: List[int]) -> List[int]:
    """
    Kasai LCP array. lcp[r] = LCP length of the suffixes sa[r] and sa[r-1];
    lcp[0] = 0. 
    
    Runs in O(N).
    """
    n = len(s)
    rank = [0] * n
    for i, p in enumerate(sa):
        rank[p] = i

    lcp = [0] * n
    h = 0
    for i in range(n):
        r = rank[i]
        if r == 0:
            h = 0
            continue
        j = sa[r - 1]
        while i + h < n and j + h < n and s[i + h] == s[j + h]:
            h += 1
        lcp[r] = h
        if h > 0:
            h -= 1
    return lcp



# multiMUM discovery


def sliding_window_min(arr: List[int], w: int) -> List[int]:
    """
    out[i] = min(arr[i : i + w]) for every valid start i, in O(len(arr)).
    """
    if w <= 0 or len(arr) < w:
        return []
    dq: collections.deque = collections.deque()
    out = [0] * (len(arr) - w + 1)
    for i, v in enumerate(arr):
        while dq and arr[dq[-1]] >= v:
            dq.pop()
        dq.append(i)
        if dq[0] <= i - w:
            dq.popleft()
        if i >= w - 1:
            out[i - w + 1] = arr[dq[0]]
    return out


def find_multimums(lanes: Sequence[str], min_len: int) -> List[Anchor]:
    """
    Find every maximal exact match that is common to all lanes and unique in each (a multiMUM), in linear time from the generalized suffix array.

    A multiMUM of length L is a block of exactly k = len(lanes) consecutive suffix-array entries such that:
        * the block covers all k lanes, once each (uniqueness + commonality),
        * L = min internal LCP of the block  (right-maximality),
        * the LCP just before and just after the block are both < L (the block is the complete occurrence set), and
        * the preceding characters are not all identical (left-maximality).
    """
    k = len(lanes)
    if k < 2 or any(len(x) == 0 for x in lanes):
        return []

    text, owner, local, alphabet_size = build_generalized(lanes)
    n = len(text)
    sa = sa_is(text, alphabet_size)
    lcp = kasai_lcp(text, sa)

    window_min = sliding_window_min(lcp, k - 1)  # min over the k-1 internal LCPs

    lane0 = lanes[0]
    anchors: List[Anchor] = []

    # Slide a width-k window over the suffix array, tracking lane coverage.
    counts = [0] * k
    distinct = 0
    separators = 0
    for j in range(k):
        o = owner[sa[j]]
        if o == -1:
            separators += 1
        else:
            if counts[o] == 0:
                distinct += 1
            counts[o] += 1

    for i in range(0, n - k + 1):
        if i > 0:
            o_out = owner[sa[i - 1]]
            if o_out == -1:
                separators -= 1
            else:
                counts[o_out] -= 1
                if counts[o_out] == 0:
                    distinct -= 1
            o_in = owner[sa[i + k - 1]]
            if o_in == -1:
                separators += 1
            else:
                if counts[o_in] == 0:
                    distinct += 1
                counts[o_in] += 1

        if separators or distinct != k:
            continue

        length = window_min[i + 1]  # min of lcp[i+1 .. i+k-1]
        if length < min_len:
            continue

        before = lcp[i]                       # lcp[0] == 0 handles i == 0
        after = lcp[i + k] if i + k < n else -1
        if before >= length or after >= length:
            continue

        positions = [0] * k
        left_maximal = False
        preceding = set()
        for m in range(i, i + k):
            g = sa[m]
            positions[owner[g]] = local[g]
            if local[g] == 0:
                left_maximal = True
            else:
                preceding.add(text[g - 1])
        if not left_maximal and len(preceding) > 1:
            left_maximal = True
        if not left_maximal:
            continue

        p0 = positions[0]
        anchors.append(Anchor(tuple(positions), length, lane0[p0:p0 + length]))

    return anchors



# Collinear chaining

def _chain_reconstruct(anchors: List[Anchor], parent: List[int], best: int) -> List[Anchor]:
    chain: List[Anchor] = []
    cur = best
    while cur != -1:
        chain.append(anchors[cur])
        cur = parent[cur]
    chain.reverse()
    return chain


def chain_pairwise(anchors: List[Anchor]) -> List[Anchor]:
    """
    Maximum-weight collinear non-overlapping chain for k == 2, in O(z log z).

    Anchor i may precede j iff it ends before j starts in both lanes:
        pos_i[d] + len_i <= pos_j[d]   for d in {0, 1}.
    Weight is anchor length. A Fenwick tree over the (compressed) lane-1 coordinate answers the prefix-maximum queries; anchors are activated in
    lane-0 end order as the sweep advances in lane-0 start order.
    """
    z = len(anchors)
    if z == 0:
        return []

    ys = sorted({a.positions[1] for a in anchors} |
                {a.positions[1] + a.length for a in anchors})
    y_rank = {y: i + 1 for i, y in enumerate(ys)}
    m = len(ys)

    tree_val = [-1] * (m + 1)
    tree_arg = [-1] * (m + 1)

    def update(pos: int, val: int, arg: int) -> None:
        while pos <= m:
            if val > tree_val[pos]:
                tree_val[pos] = val
                tree_arg[pos] = arg
            pos += pos & (-pos)

    def query(pos: int) -> Tuple[int, int]:
        best_val, best_arg = 0, -1
        while pos > 0:
            if tree_val[pos] > best_val:
                best_val = tree_val[pos]
                best_arg = tree_arg[pos]
            pos -= pos & (-pos)
        return best_val, best_arg

    by_start = sorted(range(z), key=lambda i: (anchors[i].positions[0],
                                               anchors[i].positions[1]))
    by_end = sorted(range(z), key=lambda i: anchors[i].positions[0] + anchors[i].length)

    dp = [0] * z
    parent = [-1] * z
    best_idx = by_start[0]
    best_total = 0
    ptr = 0

    for idx in by_start:
        a = anchors[idx]
        start0 = a.positions[0]
        start1 = a.positions[1]
        while ptr < z:
            bi = by_end[ptr]
            b = anchors[bi]
            if b.positions[0] + b.length > start0:
                break
            update(y_rank[b.positions[1] + b.length], dp[bi], bi)
            ptr += 1
        prev_val, prev_arg = query(y_rank[start1])
        dp[idx] = prev_val + a.length
        parent[idx] = prev_arg
        if dp[idx] > best_total:
            best_total = dp[idx]
            best_idx = idx

    return _chain_reconstruct(anchors, parent, best_idx)


def chain_multi(anchors: List[Anchor]) -> List[Anchor]:
    """
    Maximum-weight collinear non-overlapping chain for any number of lanes, by exact dynamic programming. Anchors are multiMUMs (unique per lane)
    """
    z = len(anchors)
    if z == 0:
        return []

    order = sorted(anchors, key=lambda a: (a.positions[0], a.positions))
    dim = len(order[0].positions)
    dp = [0] * z
    parent = [-1] * z
    best_idx = 0
    best_total = 0

    for j in range(z):
        aj = order[j]
        best_prev = 0
        best_prev_i = -1
        jp = aj.positions
        for i in range(j):
            ai = order[i]
            ip = ai.positions
            ilen = ai.length
            if all(ip[d] + ilen <= jp[d] for d in range(dim)):
                if dp[i] > best_prev:
                    best_prev = dp[i]
                    best_prev_i = i
        dp[j] = best_prev + aj.length
        parent[j] = best_prev_i
        if dp[j] > best_total:
            best_total = dp[j]
            best_idx = j

    return _chain_reconstruct(order, parent, best_idx)


def select_anchor_chain(anchors: List[Anchor]) -> List[Anchor]:
    if not anchors:
        return []
    if len(anchors[0].positions) == 2:
        return chain_pairwise(anchors)
    return chain_multi(anchors)



# ED string construction

def make_ed_segment(variants: Sequence[str]) -> str:
    """
    Render one ED segment from the per-lane substrings of a gap.

    Equal variants collapse to a deterministic segment (empty -> empty string).
    Otherwise a nondeterministic {...} group is emitted, deduplicated and order-preserving, with the empty word rendered as an empty variant.
    """
    unique: List[str] = []
    seen = set()
    for v in variants:
        if v not in seen:
            seen.add(v)
            unique.append(v)
    if len(unique) == 1:
        return unique[0]
    return "{" + ",".join(unique) + "}"


def build_ed(
    lanes: Sequence[str],
    min_anchor_len: int,
    max_depth: int,
    stats: Dict[str, int],
    depth: int = 0,
) -> str:
    """
    Recursively decompose k aligned lanes into an ED string.

    Anchors common to all lanes become deterministic segments; the gaps between them are decomposed the same way one level deeper, so matches that are only
    locally unique still become anchors. A gap with no common anchor (or an empty lane, or exhausted depth) collapses to one nondeterministic segment.
    """
    k = len(lanes)
    if k == 0:
        return ""
    if k == 1:
        return lanes[0]

    first = lanes[0]
    if all(lane == first for lane in lanes):
        return first

    longest = max(len(lane) for lane in lanes)
    if (depth >= max_depth
            or longest < min_anchor_len
            or any(len(lane) == 0 for lane in lanes)):
        return make_ed_segment(lanes)

    anchors = find_multimums(lanes, min_anchor_len)
    if depth == 0:
        stats["candidates"] = len(anchors)

    chain = select_anchor_chain(anchors)
    if not chain:
        return make_ed_segment(lanes)

    parts: List[str] = []
    cursors = [0] * k
    for anchor in chain:
        gap = [lanes[i][cursors[i]:anchor.positions[i]] for i in range(k)]
        segment = build_ed(gap, min_anchor_len, max_depth, stats, depth + 1)
        if segment:
            parts.append(segment)

        parts.append(anchor.text)
        stats["chain_anchors"] += 1
        stats["total_anchor_length"] += anchor.length

        for i in range(k):
            cursors[i] = anchor.positions[i] + anchor.length

    tail = [lanes[i][cursors[i]:] for i in range(k)]
    segment = build_ed(tail, min_anchor_len, max_depth, stats, depth + 1)
    if segment:
        parts.append(segment)

    return "".join(parts)



# Per-block processing

def deduplicate_sequences(records: Sequence[FastaRecord]) -> List[str]:
    """Distinct sequences, first-occurrence order preserved."""
    seen = set()
    seqs: List[str] = []
    for r in records:
        if r.seq not in seen:
            seen.add(r.seq)
            seqs.append(r.seq)
    return seqs


def process_block(
    block_id: str,
    records: Sequence[FastaRecord],
    min_anchor_len: int,
    max_depth: int,
) -> Tuple[str, str, Dict[str, int]]:
    """Process one synteny block. Returns (header, ed_string, stats)."""
    seqs = deduplicate_sequences(records)

    stats = {
        "input_sequences": len(records),
        "unique_sequences": len(seqs),
        "candidates": 0,
        "chain_anchors": 0,
        "total_anchor_length": 0,
    }

    if not seqs:
        return block_id, "", stats

    if len(seqs) == 1:
        stats["chain_anchors"] = 1 if seqs[0] else 0
        stats["total_anchor_length"] = len(seqs[0])
        return block_id, seqs[0], stats

    ed = build_ed(seqs, min_anchor_len, max_depth, stats)

    return block_id, ed, stats



# CLI

def group_by_block(records: Sequence[FastaRecord]) -> Dict[str, List[FastaRecord]]:
    groups: Dict[str, List[FastaRecord]] = collections.defaultdict(list)
    for r in records:
        groups[r.block_id].append(r)
    return dict(groups)


def _format_stats_line(block_id: str, stats: Dict[str, int]) -> str:
    return (
        f"[block {block_id}] "
        f"seqs={stats['input_sequences']} "
        f"unique={stats['unique_sequences']} "
        f"candidates={stats['candidates']} "
        f"anchors={stats['chain_anchors']} "
        f"anchor_len={stats['total_anchor_length']}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Construct ED strings from FASTA sequences grouped by synteny-block identifier."
    )
    parser.add_argument("-i", "--input", help="Input FASTA file.")
    parser.add_argument("-o", "--output",
                        help="Output FASTA file (one ED string per block).")
    parser.add_argument(
        "--min-anchor-len", type=int, default=8,
        help="Minimum length of a common exact anchor. Smaller values yield finer, more deterministic ED strings but run slower. Default: 8.",
    )
    parser.add_argument(
        "--max-depth", type=int, default=6,
        help="Maximum recursive gap-refinement depth. Default: 6.",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1,
        help="Number of parallel worker processes. Default: 1.",
    )
    parser.add_argument("--stats", default=None,
                        help="Optional path for per-block statistics (TSV).")

    parser.add_argument("--max-candidates-per-window", type=int, default=None,
                        help=argparse.SUPPRESS)
    parser.add_argument("--max-total-candidates", type=int, default=None,
                        help=argparse.SUPPRESS)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.input or not args.output:
        parser.error("the following arguments are required: -i/--input, -o/--output")
    if args.min_anchor_len < 1:
        parser.error("--min-anchor-len must be at least 1")
    if args.max_depth < 0:
        parser.error("--max-depth must be non-negative")

    records = read_fasta(args.input)
    if not records:
        print("No FASTA records found.", file=sys.stderr)
        sys.exit(1)

    groups = group_by_block(records)
    print(f"Read {len(records)} sequences.", file=sys.stderr)
    print(f"Found {len(groups)} synteny blocks.", file=sys.stderr)

    items = sorted(groups.items(), key=lambda kv: kv[0])
    results: List[Tuple[str, str]] = []
    stats_rows: List[Tuple[str, Dict[str, int]]] = []

    if args.threads <= 1:
        for block_id, block_records in items:
            header, ed, stats = process_block(
                block_id, block_records, args.min_anchor_len,
                args.max_depth,
            )
            results.append((header, ed))
            stats_rows.append((block_id, stats))
            print(_format_stats_line(block_id, stats), file=sys.stderr)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as ex:
            future_to_block = {
                ex.submit(process_block, block_id, block_records,
                          args.min_anchor_len, args.max_depth): block_id
                for block_id, block_records in items
            }
            collected: Dict[str, Tuple[str, str, Dict[str, int]]] = {}
            for fut in concurrent.futures.as_completed(future_to_block):
                block_id = future_to_block[fut]
                header, ed, stats = fut.result()
                collected[block_id] = (header, ed, stats)
                print(_format_stats_line(block_id, stats), file=sys.stderr)
            for block_id, _ in items:
                header, ed, stats = collected[block_id]
                results.append((header, ed))
                stats_rows.append((block_id, stats))

    write_fasta(results, args.output)

    if args.stats is not None:
        with open(args.stats, "w", encoding="utf-8") as f:
            f.write("block_id\tinput_sequences\tunique_sequences\tcandidates\t"
                    "chain_anchors\ttotal_anchor_length\n")
            for block_id, stats in stats_rows:
                f.write(f"{block_id}\t{stats['input_sequences']}\t"
                        f"{stats['unique_sequences']}\t{stats['candidates']}\t"
                        f"{stats['chain_anchors']}\t{stats['total_anchor_length']}\n")

    print(f"Wrote ED strings to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
