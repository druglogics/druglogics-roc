# druglogics-roc-generator

A Shiny R document to generate ROC curves based on Drabme ensemble-wise synergy 
output results. You can run the example using Rstudio or from the command line:
```
Rscript -e "library(rmarkdown); rmarkdown::run(file="roc_gen.Rmd", shiny_args = list(port = 3333))"
```