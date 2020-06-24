rm(list=ls())

# Load packages
require(knitr)
require(markdown)

setwd("S:\\Documents\\AuxMixtures\\Demo\\")

knit("SimSSVS.Rmd")
markdownToHTML("SimSSVS.md","SimSSVS.html")

knit("SimFA.Rmd")
markdownToHTML("SimFA.md","SimFA.html")

knit("SimPCR.Rmd")
markdownToHTML("SimPCR.md","SimPCR.html")
