# FireSpread
Simulate fire spread with a cellular automata approach.

## Building and installing
For building:
```sh
Rscript -e "Rcpp::compileAttributes()"
R CMD build .
```
That should generate a file called `FireSpread_1.1.tar.gz`.

For installing:
```sh
R CMD INSTALL FireSpread_1.0.tar.gz
```
