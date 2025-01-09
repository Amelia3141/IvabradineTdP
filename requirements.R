install.packages(c(
  "shiny",
  "reticulate",
  "shinythemes",
  "DT",
  "plotly"
))

# Install Python dependencies through reticulate
reticulate::py_install(c(
  "torch",
  "transformers",
  "sentence-transformers",
  "biopython",
  "flask"
))
