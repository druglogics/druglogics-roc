# druglogics-roc-generator

An R markdown document with Shiny components to generate ROC curves based on 
Drabme ensemble-wise synergy output results. You can run the example using 
Rstudio or from the command line:
```
Rscript -e 'library(rmarkdown); rmarkdown::run(file = paste0(getwd(), "/roc_gen.Rmd"), shiny_args = list(host = "127.0.0.1", port = 3333, launch.browser = TRUE))'
```