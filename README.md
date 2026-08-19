# Ja-be-Ja Graph Partitioning

Java implementation of the **Ja-be-Ja** distributed graph partitioning algorithm, extended
with a simulated-annealing acceptance criterion. Originally a lab assignment for the ID2222
Data Mining course at KTH (Group 4).

## What it does

Ja-be-Ja partitions a graph into `k` balanced partitions while minimizing the edge cut
(edges between nodes in different partitions), without any central coordinator: each node
repeatedly samples a set of candidate partners — its direct neighbors and/or a random sample
of the whole graph — and swaps its partition ("color") with the partner if doing so improves
a local energy function. A temperature parameter `T` controls how eagerly bad swaps are
accepted early on and cools down over the rounds, similar to simulated annealing.

The project explores three variants of the acceptance rule (implemented in
[`Jabeja.java`](src/main/java/se/kth/jabeja/Jabeja.java)):

- **Task 1** — the original Ja-be-Ja deterministic acceptance rule.
- **Task 2** — a simulated-annealing variant with periodic temperature restarts, tuned across
  several cooling schedules (`try2`, `try3`).
- **Optional task** — a probabilistic (Metropolis-style) acceptance rule.

## Project structure

- [`src/main/java/se/kth/jabeja/`](src/main/java/se/kth/jabeja) — algorithm implementation
  (`Jabeja.java`), graph model (`Node.java`), CLI/config handling, and graph I/O.
- [`graphs/`](graphs) — input graphs used for the experiments (`3elt`, `add20`, `facebook`,
  `twitter`, synthetic graphs, ...).
- `task1/`, `task2/`, `task2_optional/` — result logs (edge-cut, swaps, migrations per round)
  produced by each experiment.
- `graph_*.png` — plotted results per graph/task.
- [`compile.sh`](compile.sh), [`run.sh`](run.sh), [`plot.sh`](plot.sh) — build, run, and
  plotting scripts.

## Build & run

```bash
./compile.sh   # mvn clean install
./run.sh -graph ./graphs/3elt.graph -rounds 1000 -nodeSelectionPolicy HYBRID
```

Available CLI options (see [`CLI.java`](src/main/java/se/kth/jabeja/io/CLI.java)):

| Option | Description | Default |
|---|---|---|
| `-graph` | Input graph file | `./graphs/ws-250.graph` |
| `-rounds` | Number of rounds | `1000` |
| `-numPartitions` | Number of partitions | `4` |
| `-nodeSelectionPolicy` | `RANDOM`, `LOCAL`, or `HYBRID` | `HYBRID` |
| `-graphInitColorSelectionPolicy` | `RANDOM`, `ROUND_ROBIN`, or `BATCH` | `ROUND_ROBIN` |
| `-temp` | Initial simulated-annealing temperature | `2` |
| `-delta` | Cooling rate | `0.003` |
| `-alpha` | Energy function exponent | `2` |
| `-randNeighborsSampleSize` | Neighbors sampled per round | `3` |
| `-uniformRandSampleSize` | Random graph sample size per round | `6` |

Results are written under `./output` as tab-separated round / edge-cut / swaps / migrations
logs. Plot a result file with:

```bash
./plot.sh output/<result-file>.txt
```

## Requirements

- Java 7+
- Maven
- gnuplot (for plotting)
