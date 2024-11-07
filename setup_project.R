# Crear estructura de directorios
create_project_structure <- function() {
  dirs <- c(
    "scripts/preprocessing",
    "scripts/imputation",
    "scripts/analysis",
    "scripts/visualization",
    "scripts/utils",
    "data/raw",
    "data/processed",
    "results/plots",
    "results/tables",
    "docs"
  )
  
  # Crear directorios
  sapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  
  # Crear archivos placeholder
  placeholders <- c(
    "data/raw/.gitkeep",
    "data/processed/.gitkeep",
    "results/plots/.gitkeep",
    "results/tables/.gitkeep"
  )
  
  sapply(placeholders, file.create)
  
  cat("Estructura de directorios creada con éxito.\n")
}

# Crear .gitignore
create_gitignore <- function() {
  gitignore_content <- c(
    "# R files",
    ".Rproj.user",
    ".Rhistory",
    ".RData",
    ".Ruserdata",
    "",
    "# Data",
    "/data/raw/*",
    "/data/processed/*",
    "!/data/raw/.gitkeep",
    "!/data/processed/.gitkeep",
    "",
    "# Results",
    "/results/*",
    "!/results/plots/.gitkeep",
    "!/results/tables/.gitkeep",
    "",
    "# System files",
    ".DS_Store",
    "Thumbs.db",
    "",
    "# Log files",
    "*.log"
  )
  
  writeLines(gitignore_content, ".gitignore")
  cat("Archivo .gitignore creado con éxito.\n")
}

# Función principal de configuración
setup_project <- function() {
  cat("Iniciando configuración del proyecto...\n")
  create_project_structure()
  create_gitignore()
  cat("Configuración del proyecto completada.\n")
}
