# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install required R packages
if (!require("shiny")) install.packages("shiny")
if (!require("shinythemes")) install.packages("shinythemes")
if (!require("DT")) install.packages("DT")
if (!require("plotly")) install.packages("plotly")
if (!require("shinyjs")) install.packages("shinyjs")
if (!require("reticulate")) install.packages("reticulate")
if (!require("markdown")) install.packages("markdown")

# Configure reticulate and install Python dependencies
library(reticulate)
py_install(c("biopython", "pandas", "requests"), pip=TRUE)
