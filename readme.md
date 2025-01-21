# Bioinformatics

## Find String:

* Knuth Morris Pratt Algorithm
    * with Prefix-Function
    * with Z-Function
* Boyer Moore Algorithm
* Silly Search Algorithm

<details>
  <summary>Benchmark</summary>

```
2025-01-21T00:07:17+01:00
Running ./benchmark_bioinformatics
Run on (22 X 4500 MHz CPU s)
Load Average: 0.63, 1.07, 0.99
***WARNING*** CPU scaling is enabled, the benchmarks real time measurements may be noisy and will incur extra overhead.
***WARNING*** Library was built as DEBUG. Timings may be affected.
```

| **Benchmark**                              | **Time (ns)** |  **CPU (ns)** | **Iterations** |
|--------------------------------------------|--------------:|--------------:|---------------:|
| BM_KnuthMorrisPratt_Prefix/1000/10         |        13,714 |        13,685 |         52,155 |
| BM_KnuthMorrisPratt_Prefix/1000000/100     |    12,697,141 |    12,670,310 |             54 |
| BM_KnuthMorrisPratt_Prefix/1000000000/1000 |    1.6009e+10 |    1.4721e+10 |              1 |
| BM_KnuthMorrisPratt_Z/1000/10              |        19,633 |        19,590 |         34,159 |
| BM_KnuthMorrisPratt_Z/1000000/100          |    20,460,156 |    20,416,654 |             36 |
| BM_KnuthMorrisPratt_Z/1000000000/1000      |    2.3067e+10 |    2.2911e+10 |              1 |
| BM_BoyerMoore/1000/10                      |        27,966 |        27,866 |        120,377 |
| BM_BoyerMoore/1000000/100                  |     7,123,170 |     7,104,603 |            126 |
| BM_BoyerMoore/1000000000/1000              | 6,394,478,659 | 6,334,507,181 |              1 |
| BM_Search/1000/10                          |         9,997 |         9,974 |         69,408 |
| BM_Search/1000000/100                      |     8,431,268 |     8,408,135 |             76 |
| BM_Search/1000000000/1000                  | 9,446,965,325 | 9,404,504,676 |              1 |
| BM_StdFind/1000/10                         |         1,782 |         1,775 |        404,930 |
| BM_StdFind/1000000/100                     |     1,637,439 |     1,611,711 |            438 |
| BM_StdFind/1000000000/1000                 | 1,581,023,622 | 1,577,786,916 |              1 |

</details>

## Alignment String:

* Needleman Wunsch Algorithm
* Smith Waterman Algorithm
* Repeated Local Alignment
* Wagner Fischer Algorithm
*

<details>
  <summary>Scores Example</summary>

## Needleman-Wunsch Matrix

|   |     |   H |   E |   A |   G |   A |   W |   G |   H |   E |   E |
|---|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|
|   |   0 |  -8 | -16 | -24 | -32 | -40 | -48 | -56 | -64 | -72 | -80 |
| P |  -8 |  -2 |  -9 | -17 | -25 | -33 | -41 | -49 | -57 | -65 | -73 |
| A | -16 | -10 |  -3 |  -4 | -12 | -20 | -28 | -36 | -44 | -52 | -60 |
| W | -24 | -18 | -11 |  -6 |  -7 | -15 |  -5 | -13 | -21 | -29 | -37 |
| H | -32 | -14 | -18 | -13 |  -8 |  -9 | -13 |  -7 |  -3 | -11 | -19 |
| E | -40 | -22 |  -8 | -16 | -16 |  -9 | -12 | -15 |  -7 |   3 |  -5 |
| A | -48 | -30 | -16 |  -3 | -11 | -11 | -12 | -12 | -15 |  -5 |   2 |
| E | -56 | -38 | -24 | -11 |  -6 | -12 | -14 | -15 | -12 |  -9 |   1 |

## Smith-Waterman Matrix

|   |   |  H |  E |  A |  G |  A |  W |  G |  H |  E |  E |
|---|--:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|   | 0 |    |    |    |    |    |    |    |    |    |    |
| P |   |    |    |    |    |    |    |    |    |    |    |
| A |   |    |    |  5 |    |  5 |    |    |    |    |    |
| W |   |    |    |    |  2 |    | 20 | 12 |  4 |    |    |
| H |   | 10 |  2 |    |    |    | 12 | 18 | 22 | 14 |  6 |
| E |   |  2 | 16 |  8 |    |    |  4 | 10 | 18 | 28 | 20 |
| A |   |    |  8 | 21 | 13 |  5 |    |  4 | 10 | 20 | 27 |
| E |   |    |  6 | 13 | 18 | 12 |  4 |    |  4 | 16 | 26 |

## Repeated Local Alignment Matrix

|   |   |  H |  E |  A |  G |  A |  W |  G |  H |  E |  E |
|---|--:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|   | 0 |    |    |    |  1 |  1 |  1 |  1 |  1 |  3 |  9 |
| P |   |    |    |    |  1 |  1 |  1 |  1 |  1 |  3 |  9 |
| A |   |    |    |  5 |  1 |  6 |  1 |  1 |  1 |  3 |  9 |
| W |   |    |    |    |  2 |  1 | 21 | 13 |  5 |  3 |  9 |
| H |   | 10 |  2 |    |  1 |  1 | 13 | 19 | 23 | 15 |  9 |
| E |   |  2 | 16 |  8 |  1 |  1 |  5 | 11 | 19 | 29 | 21 |
| A |   |    |  8 | 21 | 13 |  6 |  1 |  5 | 11 | 21 | 28 |
| E |   |    |  6 | 13 | 18 | 12 |  4 |  1 |  5 | 17 | 27 |

## Wagner-Fischer Matrix

|   |   | H | E | A | G | A | W | G | H | E |  E |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|---:|
|   | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
| P | 1 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
| A | 2 | 2 | 2 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |  9 |
| W | 3 | 3 | 3 | 3 | 3 | 4 | 4 | 5 | 6 | 7 |  8 |
| H | 4 | 3 | 4 | 4 | 4 | 4 | 5 | 5 | 5 | 6 |  7 |
| E | 5 | 4 | 3 | 4 | 5 | 5 | 5 | 6 | 6 | 5 |  6 |
| A | 6 | 5 | 4 | 3 | 4 | 5 | 6 | 6 | 7 | 6 |  6 |
| E | 7 | 6 | 5 | 4 | 4 | 5 | 6 | 7 | 7 | 7 |  6 |

</details>