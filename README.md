Análisis de Alelos HLA en COVID-19
Descripción
Este proyecto analiza la relación entre alelos HLA y la severidad de COVID-19, incluyendo análisis de hospitalización y casos asintomáticos.
Estructura del Proyecto
hla-covid-analisis/
├── scripts/              # Scripts de análisis
│   ├── preprocessing/    # Preprocesamiento de datos
│   ├── imputation/      # Imputación HLA
│   ├── analysis/        # Análisis estadístico
│   ├── visualization/   # Visualización
│   └── utils/          # Utilidades
├── data/                # Datos
│   ├── raw/            # Datos sin procesar
│   └── processed/      # Datos procesados
├── results/            # Resultados
├── plots/              # Visualizaciones
└── docs/              # Documentación
Requisitos

R versión 4.0.0 o superior
Paquetes R requeridos:

HIBAG
data.table
ggplot2
[otros paquetes...]



Instalación

Clonar el repositorio:

git clone https://github.com/davidgarvalfer/hla-covid-analisis.git
cd hla-covid-analisis

Instalar dependencias:

source("scripts/utils/config.R")
Uso
source("scripts/main.R")
Configuración
Modificar el archivo config.yml según sea necesario:

Rutas de datos
Parámetros de análisis
Configuración de resultados

Resultados
Los resultados se guardan en:

results/: Archivos de resultados estadísticos
plots/: Visualizaciones generadas

Autor
David García Valentín-Fernández
