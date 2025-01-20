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

