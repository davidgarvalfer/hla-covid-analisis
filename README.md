Análisis HLA-COVID
Análisis estadístico de la relación entre alelos HLA y la severidad de COVID-19.
Estructura del Proyecto
Copyhla-covid-analisis/
├── scripts/
│   ├── preprocessing/      # Preprocesamiento de datos
│   ├── imputation/        # Imputación HLA
│   ├── analysis/          # Análisis estadístico
│   ├── visualization/     # Visualización de resultados
│   └── utils/             # Utilidades y configuración
├── data/
│   ├── raw/              # Datos originales
│   └── processed/        # Datos procesados
├── results/
│   ├── plots/            # Gráficos generados
│   └── tables/           # Tablas de resultados
└── docs/                 # Documentación adicional
Requisitos

R versión 4.0.0 o superior
Paquetes R necesarios:

HIBAG
data.table
ggplot2
corrplot
pROC
[etc...]



Instalación

Clonar el repositorio
Ejecutar source("scripts/utils/config.R") para instalar dependencias
Ejecutar source("scripts/main.R") para realizar el análisis completo

Uso
rCopy# Ejecutar análisis completo
source("scripts/main.R")

# O ejecutar pasos individuales
source("scripts/preprocessing/data_preparation.R")
source("scripts/imputation/imputation_functions.R")
source("scripts/analysis/statistical_analysis.R")
source("scripts/visualization/plotting_functions.R")
Resultados
Los resultados se guardan en:

Gráficos: results/plots/
Tablas: results/tables/

Contacto
David García