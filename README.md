# druglogics-roc

A [Shiny](https://shiny.rstudio.com/) app that generates ROC curves based on Drabme's [ensemble-wise synergy output](https://druglogics.github.io/druglogics-doc/drabme-install.html#drabme-output) results.

## Install

The following libraries must be installed before using the app:

- [shiny](https://cran.r-project.org/web/packages/shiny/index.html)
- [emba](https://cran.r-project.org/web/packages/emba/index.html)
- [usefun](https://cran.r-project.org/web/packages/usefun/index.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
- [DT](https://cran.r-project.org/web/packages/DT/index.html)
- [plotly](https://github.com/ropensci/plotly)

## Run

You can test the app without cloning the repo: `shiny::runGitHub("druglogics/druglogics-roc")`.

Otherwise, clone the repo and run: `shiny::runApp()`

## Example Input

The app generates ROC curves using the [example test files](https://github.com/druglogics/druglogics-roc/tree/master/examples).
The user can also upload his own **input files** to generate a ROC curve.
The input files must be:

- An *observed synergies* file to test against the synergies found from `Drabme`'s simulations.
These are the *true positive* synergies.
See [example file](https://github.com/druglogics/druglogics-roc/blob/master/examples/observed_synergies).
- An [ensemble-wise synergy scores](https://druglogics.github.io/druglogics-doc/drabme-install.html#drabme-output) result file (most common use case is to run the Gitsbe and then Drabme using the `druglogics-synergy` module to get this output file).
See respective [example file](https://github.com/druglogics/druglogics-roc/blob/master/examples/ensemble_synergies).

In the case that two ensemble-wise synergy score files are used (usually the second represents the results of a Gitsbe analysis where the models are [trained to a proliferation profile](https://druglogics.github.io/druglogics-doc/training-data.html#unperturbed-condition---globaloutput-response) - i.e. *random* models) we can merge the synergy scores to a **combined-classifier** output and produce thus a new combined ROC curve.

The app's interface allows the user to choose which **output** result will be used to generate a ROC curve.
The options for the single output results are `original` and `random` and for the combined classifier output `original-random` and `original+random`.
