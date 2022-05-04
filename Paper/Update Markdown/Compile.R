library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))
getwd()

bookdown::render_book('index.Rmd', "bookdown::pdf_document2", new_session = T)

