## 2026-11-01 Estructura general

**Corpus.** BIOFACQUIM V2
Descargar la base de datos BioFacQuim desde: (https://doi.org/10.6084/m9.figshare.11312702), son compuestos en formato SDF con licencia CCO. Son metabolitos secundarios de origen mexicano con anotaciones taxonómicas y procedencia por compuesto.

**Criterios de clasificación** 
Realizar la deduplicación (eliminar copias redundantes de la base de datos) con el primer bloque del InChIKey de 14 caracteres (Es la información que solo incluye la conectividad y el esqueleto molecular). Se debe reportar los que tienen estereoquímica y los que solo tienen la conectividad molecular.

**Masa molecular** 
Los límites inferiores y superiores serán de 50 a 5000 Da, este último ya que existen glucósidos que alcanzas masas altas rápidamente. Se debe reportar cuántos metabolitos quedan fuera de estos rangos.

#2026-11-1 Crear el entorno y clonar el repositorio

1) Se clona y hace fork al repositorio orignal

2) Se crea y cambia el flujo de trabajo a la rama "etapa1-curacion"

3) Se intala y activa el entorno quimioinfo

4) Se verfica el entorno así como la existencia de todos los directorios correspondientes además de verificar los criterios del .gitignore

5) Se removió la carpeta .vscode del rastreo de git, además de añadir al .gitignore los direcotrios de resultados y figuras

#2026-11-02 Desarrollo y depuración del cuaderno de curación

1) Conversión y Estandarización: Se desarrolló el código para leer el archivo BIOFACQUIM_V2.sdf y convertirlo a un DataFrame de Pandas. Se implementó un pipeline de estandarización con RDKit (rdMolStandardize) para normalizar grupos, neutralizar, reionizar y seleccionar el fragmento mayor. Se generó un control visual 2D (control_estandarizacion.png).

2) Corrección lógica y aplicación de filtros: Se depuró un error lógico crítico donde el dataset se exportaba en crudo. Se creó la celda de ejecución para aplicar secuencialmente los filtros predefinidos:

Filtro de elementos orgánicos (remoción de metales/inorgánicos).

Filtro de peso molecular (50 - 5000 Da).

Deduplicación por esqueleto InChIKey (bloque 1 ignorando estereoquímica).

3) Mejora del reporte estadístico (Conteos): Se enriqueció la bitácora automatizada de filtrado (df_conteos) calculando métricas adicionales: moléculas descartadas, porcentaje de retención por paso y porcentaje de retención total acumulado.

4) Exportación unificada: Se configuró el script para crear automáticamente los subdirectorios si no existen (os.makedirs) en las rutas personalizadas (torresC1). Se integró la exportación de resultados al final del cuaderno en una sola celda para permitir su ejecución como script, generando tres archivos:

torres_conteos_BIOFACQUIM_V2.csv (Bitácora de filtrado).

torresBIOFACQUIM_V2curado.csv (Dataset tabular para ML y análisis).

torresBIOFACQUIM_V2curado.sdf (Dataset estructural conservando metadatos para visualización 3D y Docking).

##Fecha en que se descarga 2026-11-1 Versión del conjunto BIOFACQUIM_V2.sdf sin curar. Licencia CCO

#2026-11-01 | Curación del Corpus BIOFACQUIM V2
Curación de Datos
Fuente del Corpus:

Dataset: BIOFACQUIM v2 (metabolitos secundarios de origen mexicano con anotaciones taxonómicas y de procedencia geográfica).

Formato original: .sdf bajo licencia Creative Commons Zero (CC0).

Identidad digital / Enlace: DOI: 10.6084/m9.figshare.11312702.

Criterios de Deduplicación:

Extracción y agrupar por el primer bloque del InChIKey (14 caracteres), el cual codifica únicamente la conectividad y el esqueleto molecular base.

Reporte comparativo de registros con estereoquímica explícita definida frente a estructuras limitadas a conectividad molecular plana.

Filtro de Masa Molecular (MW):

Rango permitido: 50 Da – 5000 Da.

Nota de corte: El límite superior de 5000 Da responde a la presencia potencial de glicósidos y conjugados complejas de alto peso molecular. Se registrará la cantidad exacta de estructuras excluidas por debajo y por encima de estos umbrales.

## 8 de Noviembre de 2026 - Perfilado quimioinformático y análisis 

### 1. Consolidación de Datasets
- Se unificaron 3 conjuntos de datos en un solo DataFrame de trabajo:
  - Productos naturales (`curado.csv`): 503 moléculas.
  - Fármacos aprobados (`referencia_aprobados.csv`): 2,156 moléculas.
  - Moléculas aleatorias (`referencia_aleatoria.csv`): 503 moléculas.
- Se validaron y estandarizaron las columnas SMILES, procesando un total de 3,162 moléculas estructuralmente válidas bajo la variable `SMILES_UNIFICADO`.

### 2. Extracción de Características (RDKit)
- **Descriptores:** Se calcularon métricas clave (MW, cLogP, TPSA, HBD, HBA, RotB, Fracción Csp3, Anillos, Estereocentros, Átomos Pesados y QED) para el corpus total.
- **Fingerprints:** Se obtuvieron los Morgan Fingerprints utilizando un radio de 2 y un vector de 2048 bits.
- **Entregables:** Se exportó exitosamente el archivo `data/processed/torresprocessed/torresF2/descriptores.csv`, el cual incluye la matriz completa para las 503 moléculas curadas (73 columnas).

### 3. Comparación Estadística No Paramétrica
- Se evaluó el tamaño del efecto (Delta de Cliff) y se corrieron pruebas de Mann-Whitney U para determinar diferencias significativas entre el grupo de productos naturales y los grupos de control.
- Los valores $p$ fueron sometidos a corrección de Holm para comparaciones múltiples.
- **Entregables:**
  - Archivo de resultados: `results/torresresults/torresF2/comparacion_descriptores.csv` (documentando 22 comparaciones válidas).
  - Registro de exclusión: `descriptores_excluidos.csv` (para descriptores con muestra insuficiente).

### 4. Análisis de Componentes Principales (PCA)
- Se aplicó un ajuste con `StandardScaler` de forma global para los tres conjuntos, garantizando que el espacio químico de referencia no tuviera desplazamientos artificiales.
- **Resultados del PCA:**
  - **PC1 (59.8% de la varianza):** Dominado fuertemente por TPSA, MW, HBA, y conteo de Átomos Pesados.
  - **PC2 (11.2% de la varianza):** Definido principalmente por variables de rigidez estructural (Fracción Csp3) y lipofilicidad (cLogP).
