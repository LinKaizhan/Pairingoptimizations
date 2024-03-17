# Pairing optimizations for isogeny-based cryptosystems

A proof-of-concept implemtation using [SageMath](https://www.sagemath.org) of version 10.1.
This implementation contains four algorithms in isogeny-based cryptosystems:

1. Full torsion point verification in "verify_ftp.sage"
2. Full torsion point verification with known $[\lambda^{-1}]Q$ in "verify_ftpl.sage"
3. Supersingularity verification in "verify_sup.sage".
4. Torsion basis generation in "gen_basis.sage".

Each algorithm count and output the number of field operations during the computation.
To execute the algorithm, one can use the command:
```
sage "algorithm_name.sage"
```
algorithm_name can be verify_ftp, verify_ftpl, verify_sup or gen_basis.
The performance of torsion basis generation is not in constant time, here we execute 10000 benchmarks and take the average cost. It should be noticed that the file "gen_basis.sage" will run around 15 minites.