setwd("D:/Github/HPCredit/Paper/Rmarkdown")

bookdown::render_book('index.Rmd', "bookdown::pdf_document2", new_session = T)
