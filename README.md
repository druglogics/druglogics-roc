# druglogics-roc-generator

An R markdown document with an embeded Shiny app to generate ROC curves based 
on the `Drabme` ensemble-wise synergy output results. If you are using Rstudio,
try this command:  
`rmarkdown::run(file = paste0(getwd(), "/roc_gen.Rmd"), shiny_args = list(host = "127.0.0.1", port = 3333, launch.browser = TRUE))`, or from the command line:  
`Rscript -e 'library(rmarkdown); rmarkdown::run(file = paste0(getwd(), "/roc_gen.Rmd"), shiny_args = list(host = "127.0.0.1", port = 3333, launch.browser = TRUE))'`
