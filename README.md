# Construction of Elastic-Degenerate Strings from Synteny Blocks

## Overview

This program constructs one elastic-degenerate string (ED string) for each synteny block represented by a family of FASTA sequences.

FASTA records are grouped by the block identifier contained before the first `|` character in the header. Within each block, duplicate complete sequences are removed. The remaining sequences are then decomposed into shared deterministic regions and variable nondeterministic regions.

The implementation uses generalized suffix arrays, LCP arrays, common unique maximal exact matches (multiMUMs), maximum-weight anchor chaining, and recursive gap refinement.

Only the Python standard library is required.

---

## ED String Format

Deterministic regions are written as ordinary strings.

Nondeterministic regions are written as comma-separated alternatives inside braces.

Example:

```text
{A,C,}GAAT{AT,A}ATT
```

In this example:

- `{A,C,}` contains the alternatives `A`, `C`, and the empty string;
- `GAAT` and `ATT` are deterministic regions;
- `{AT,A}` contains the alternatives `AT` and `A`.

The empty string is represented by an empty field inside a nondeterministic segment.

Examples:

```text
{A,C,}
{,A,C}
{A,,C}
```

All three contain the empty-string alternative.

---

## Input Format

The input must be a FASTA file containing sequences from one or more synteny blocks.

The block identifier is extracted from the part of the FASTA header before the first `|`.

Example:

```text
>5653|GCF_000005845.2_ASM584v2_genomic.fna|NC_000913.3:1-5593(-)
ACTGACTGACTG
>5653|GCF_000008865.2_ASM886v2_genomic.fna|NC_002695.2:10-22(+)
ACTGACCGACTG
>7812|GCF_000005845.2_ASM584v2_genomic.fna|NC_000913.3:100-112(+)
GGTTAACCGGTT
```

The first two records belong to block `5653`. The third record belongs to block `7812`.

Multiline FASTA records are supported. All sequences are converted to uppercase.

---

## Output Format

The output is a FASTA-style file containing one ED string per synteny block.

Example:

```text
>5653
ACTG{A,AC}CCGACTG
>7812
GGTTAACCGGTT
```

Each ED string is written on a single line.

---

## Requirements

- Python 3.9 or newer
- No external Python packages are required

The program uses only modules from the Python standard library.

---

## Usage

Basic usage:

```bash
python3 construct_ed.py \
    --input input.fasta \
    --output output_ed.fasta \
    --min-anchor-len 3
```

Equivalent short form:

```bash
python3 construct_ed.py \
    -i input.fasta \
    -o output_ed.fasta \
    --min-anchor-len 3
```

Run with multiple worker processes:

```bash
python3 construct_ed.py \
    -i input.fasta \
    -o output_ed.fasta \
    --min-anchor-len 3 \
    --threads 4
```

Write per-block statistics:

```bash
python3 construct_ed.py \
    -i input.fasta \
    -o output_ed.fasta \
    --min-anchor-len 3 \
    --stats block_statistics.tsv
```

Disable output verification:

```bash
python3 construct_ed.py \
    -i input.fasta \
    -o output_ed.fasta \
    --min-anchor-len 3 \
    --no-verify
```

Run the internal self-tests:

```bash
python3 construct_ed.py --selftest
```

---

## CLI Arguments

| Argument | Description |
|----------|-------------|
| `-i`, `--input` | Input FASTA file containing sequences grouped by block identifier. |
| `-o`, `--output` | Output FASTA file containing one ED string per block. |
| `--min-anchor-len INT` | Minimum length of an exact common anchor. Default: `8`. |
| `--max-depth INT` | Maximum recursive gap-refinement depth. Default: `6`. |
| `-t INT`, `--threads INT` | Number of parallel worker processes. Default: `1`. |
| `--verify` | Verify that every output ED string spells all distinct input sequences. Enabled by default. |
| `--no-verify` | Disable ED-string verification. |
| `--stats PATH` | Write per-block statistics to a TSV file. |
| `--selftest` | Run the internal randomized correctness tests and exit. |
| `-h`, `--help` | Display the command-line help message. |

The statistics file contains the following columns:

```text
block_id
input_sequences
unique_sequences
candidates
chain_anchors
total_anchor_length
```

---

## Author

**Kiril Temelkov**
