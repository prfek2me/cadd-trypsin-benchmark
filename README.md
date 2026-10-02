# CADD Trypsin Benchmark: Redocking & MD Pipeline Automático

Un flujo de trabajo riguroso, reproducible y completamente automatizado (end-to-end) para validación de complejos proteína-ligando. Este proyecto resuelve la variabilidad paramétrica y la fricción manual en los estudios de Diseño de Fármacos Asistido por Computadora (CADD), implementando el sistema Benzamidina-Tripsina (PDB: 3PTB) como caso de estudio metodológico.

## Motivación

Históricamente, los pipelines de CADD y Dinámica Molecular presentan deficiencias críticas de reproducibilidad:
1. Fragmentación de los pasos estructurales (dependencia de clics manuales e interfaces gráficas rotas).
2. Colisión en las reglas de mezcla al combinar campos de fuerza heterogéneos.
3. Tratamiento ingenuo de iones estructurales clave al preparar la apoproteína.

Este repositorio encapsula una arquitectura de software *in silico* mediante **Snakemake**, desterrando el preprocesamiento ad-hoc y garantizando trazabilidad desde el archivo cristalográfico crudo hasta la evaluación termodinámica del acoplamiento.

## Arquitectura Científica y Ecosistema

Para asegurar reproducibilidad estricta y blindar el flujo contra conflictos de dependencias, el pipeline emplea una aproximación unificada sobre el estándar GROMACS/AMBER:

* **Estructura Receptora:** `PDBFixer` (retiene Ca2+ estabilizador y asigna protonación biológica a pH 7.4).
* **Topología Inyectada:** RDKit/Meeko extrae y parametriza las coordenadas cristalográficas y la carga catiónica (+1) del ligando original.
* **Docking:** GNINA (Scoring convolucional 3D `CNNscore` + auto-box cristalográfico para redocking perfecto).
* **Topología Mecano-Cuántica:** `acpype` y AmberTools asignan cargas AM1-BCC y validan topologías compatibles con GROMACS (GAFF2) dinámicamente.
* **Dinámica Molecular:** GROMACS 2024+ configurado termodinámicamente con correcciones rigurosas (`c-rescale` y `V-rescale` en baños térmicos separados).
* **Análisis In-Memory:** Procesamiento de trayectoria directo sobre binarios `.xtc` mediante `MDAnalysis`.

*(Nota sobre binarios compilados: Por estándares de control de versiones HPC, el ejecutable de GNINA no se distribuye en este repositorio y debe descargarse durante la inicialización).*

## Instalación del Entorno (Linux / WSL2)

El entorno se administra mediante Mamba/Conda. 

```bash
# 1. Clonar repositorio
git clone https://github.com/prfek2me/cadd-trypsin-benchmark.git
cd cadd-trypsin-benchmark

# 2. Instalar el ecosistema
mamba env create -f environment.yaml

# Nota: Si `mamba activate` falla, inicializa tu shell o usa conda:
# ~/miniforge3/bin/mamba shell init --shell bash && source ~/.bashrc
conda activate cadd-env

# 3. Preparar binario standalone de GNINA
mkdir -p bin
wget https://github.com/gnina/gnina/releases/download/v1.1/gnina -O bin/gnina
chmod +x bin/gnina
```