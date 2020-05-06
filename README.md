# druglogics-roc

A [Shiny](https://shiny.rstudio.com/) app that generates ROC and PR curves based on Drabme's [ensemble-wise synergy output](https://druglogics.github.io/druglogics-doc/drabme-install.html#drabme-output) results.

## Install

The following libraries must be installed before using the app:

- [shiny](https://cran.r-project.org/web/packages/shiny/index.html)
- [emba](https://cran.r-project.org/web/packages/emba/index.html)
- [usefun](https://cran.r-project.org/web/packages/usefun/index.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
- [DT](https://cran.r-project.org/web/packages/DT/index.html)
- [plotly](https://github.com/ropensci/plotly)
- [PRROC](https://cran.r-project.org/web/packages/PRROC/index.html)

## Run

You can test the app without cloning the repo: `shiny::runGitHub("druglogics/druglogics-roc")`.

Otherwise, clone the repo and run: `shiny::runApp()`

## Example Input

The app generates by default a ROC (Receiver Operator Characteristic) and a PR (Precision-Recall) curve using the [example test files](https://github.com/druglogics/druglogics-roc/tree/master/examples).
The user can also upload his own **input files** to generate a new ROC&PR curve.
The input files must be:

- An *observed synergies* file to test against the synergies found from `Drabme`'s simulations.
These are the *true positive* synergies.
See [example file](https://github.com/druglogics/druglogics-roc/blob/master/examples/observed_synergies).
- An [ensemble-wise synergy scores](https://druglogics.github.io/druglogics-doc/drabme-install.html#drabme-output) result file (most common use case is to run the Gitsbe and then Drabme using the `druglogics-synergy` module to get this output file).
See respective [example file](https://github.com/druglogics/druglogics-roc/blob/master/examples/ensemble_synergies).

In the case that two ensemble-wise synergy score files are used (usually the second represents the results of a Gitsbe analysis where the models are [trained to a proliferation profile](https://druglogics.github.io/druglogics-doc/training-data.html#unperturbed-condition---globaloutput-response)) we can merge the synergy scores to a **combined-classifier** output and produce thus new combined ROC and PR curves.

The app's interface allows the user to choose which **output** result will be used to generate the curves.
The options for the single output results are `predictor1` and `predictor2` and for the combined classifier output `predictor1-predictor2` and `predictor1+predictor2`.
