<!-- README.md is generated from README.Rmd. Please edit that file -->

# pyLLS
### a Python library for missing gene value imputation using local least square algorithm

The Local Least Square (LLS) algorithm is an algorithm that is particularly effective at imputing missing values.<br>
We developed pyLLS by implementing the LLS into python framework.<br>
Our pyLLS offers more options and is significantly faster than LLS in R.<br>
pyLLS is subjected to the MIT License, which guarantees free and commercial use for industry as well as the research community.<br>

Please report bugs to Sejin Oh, at <agicic[at]naver.com> or at
<https://github.com/osj118/pyLLS/issues>.

## Installation

You can install the pyLLS from
[pypi](https://pypi.org/project/pyLLS/) with:

``` r
pip install pyLLS
```

## Example

This is basic example:
If you run 'impute_missing_gene()' without any data,<br>
it will return its description.

``` r
import pyLLS
pyLLS.impute_missing_gene()
>>> pyLLS.impute_missing_gene()

            This function estimates missing values of the specified target probes.
            # parameters
            ref (pd.DataFrame): reference data. gene x sample (n x p) DataFrame.
            target (pd.DataFrame) : Target table containing missing values. gene x sample (i x k) DataFrame.
            metric (str) : ['correlation'(default),'L1','L2']
                           Similarity metric to prioritize the probes for each target.
            maxK : maximum number of probe to be used in missing value estimation.
            useKneedle : It determines whether Kneedle algorithm should be used (True) or not (False).
                         If useKneedle==False, then maxK probes will be used to estimate missing values.
            verbose : If True, progress is reported. Otherwise, no progress is reported.
            n_jobs : Use all threads ('all') or speicified number of threads (int)
            addK = Intenger that added to Kneedle's K to prevent underfitting.
                   This will use K+addK probes to estimate missing values of a gene.
            return_probes = if true, 'target-table and mgcp' will be returned else 'target' will be returned.
            # Return
            * target : table with estimated values of missing genes that are not present in original target table.
            matrix shape will be (n x k).
            * mgcp : missing gene correlative probes. If useKneedle == True, mgcp will have R2-square column.
            # tutorial
            <------omit-------->
```

## Parameters

to be described

``` r
to be described
```

<img src="man/figures/README-parameters-1.png" width="100%" />

``` r
to be described
```

<img src="man/figures/README-parameters-2.png" width="100%" />

## Sensitivity

to be described

``` r
to be described
```

<img src="man/figures/README-sensitivity-1.png" width="100%" />
