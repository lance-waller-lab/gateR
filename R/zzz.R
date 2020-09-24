.onAttach <- function(...) {
  packageStartupMessage(paste("\nWelcome to {gateR} version ", utils::packageDescription("gateR")$Version, "\n> help(\"gateR\") # for documentation\n> citation(\"gateR\") # for how to cite\n", sep = ""), appendLF = TRUE)
}
