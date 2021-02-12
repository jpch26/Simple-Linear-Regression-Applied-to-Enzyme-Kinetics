# Erase all outputs from analysis and data simulation ---------------------

# Run this code if you want to erase all result files from analysis and data
# simulation. 

# Presentation documents and session info file names
session_info <- c("Session_Info.txt")

# Data outputs file names from analysis scripts
analysis_outputs <- list.files("data")

# Graph outputs file names from analysis scripts
graph_outputs <- list.files("graphs")

# Delete the specified files 
file.remove(session_info)
file.remove(paste0("data/", analysis_outputs))
file.remove(paste0("graphs/", graph_outputs))
