


rmarkdown::render(input = "paper/fert_fore.Rmd",
                  output_format = bookdown::pdf_document2(
                    toc=FALSE,
                    fig_caption=T,
                    keep_tex=T
                  ),
                  output_dir="paper"
                  )
