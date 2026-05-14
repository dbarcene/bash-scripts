#!/bin/bash
# File              : single_run.sh
# Author            : David Barcene <dbarcene@indicasat.org.pa>
# Date              : 30.04.2026
# Last Modified Date: 14.05.2026
# Last Modified By  : David Barcene <dbarcene@indicasat.org.pa>
# =============================================================================
# Uso     : Editar las variables de la Sección 1 y lanzar con:
#           sbatch single_run.sh
# =============================================================================
#SBATCH --job-name=freesurfer_single
#SBATCH --output=/data/username/FS_LOGS/fs_%A_%a.out
#SBATCH --error=/data/username/FS_LOGS/fs_%A_%a.err
#SBATCH --nodes=1
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=20

# =============================================================================
# SECCIÓN 1: CONFIGURACIÓN DEL ENTORNO
# =============================================================================
module purge
module load freesurfer/8.1.0-1

# =============================================================================
# SECCIÓN 2: VARIABLES DE ENTRADA — EDITAR ANTES DE LANZAR
# =============================================================================
ID_PACIENTE="paciente006"
VISITA="VISITA_1"
ID_CASO="03-006"

# Directorios base
export SUBJECTS_DIR="/data/username/RESULTADOS"
export FS_LICENSE="${FREESURFER_HOME}/license.txt"
DICOM_BASE="/data/username/${VISITA}/${ID_CASO}/DICOM/00000001"

# Subdirectorios de secuencias (ajustar si dcmunpack muestra índices distintos)
T1_DIR="${DICOM_BASE}/00000002"      # SAG3DMPRAGE
FLAIR_DIR="${DICOM_BASE}/00000003"   # SAG3DFLAIR

# Verificar que el módulo haya definido FREESURFER_HOME
if [ -z "${FREESURFER_HOME}" ]; then
    echo "ERROR: FREESURFER_HOME no está definido tras module load. Abortando."
    exit 1
fi

# Inicializar el entorno completo de FreeSurfer (necesario para segmentHA_T1.sh,
# run_samseg y demás utilidades)
source "${FREESURFER_HOME}/SetUpFreeSurfer.sh"

# Verificar que la licencia existe
if [ ! -f "${FS_LICENSE}" ]; then
    echo "ERROR: No se encontró la licencia en ${FS_LICENSE}. Abortando."
    exit 1
fi

echo "======================================================="
echo " Iniciando pipeline FreeSurfer para: ${ID_PACIENTE}"
echo " SUBJECTS_DIR : ${SUBJECTS_DIR}"
echo " FS_LICENSE   : ${FS_LICENSE}"
echo " Fecha/Hora   : $(date)"
echo "======================================================="

# =============================================================================
# PASO 1: VOLUMETRÍA GLOBAL — recon-all (T1)
# =============================================================================
echo ""
echo "--- PASO 1: recon-all ---"
echo "Inicio: $(date)"

recon-all -all \
    -i "${T1_DIR}/00000001" \
    -s "${ID_PACIENTE}" \
    -sd "${SUBJECTS_DIR}" \
    -threads 20

echo "recon-all finalizado: $(date)"

# =============================================================================
# PASO 2: SEGMENTACIÓN DEL HIPOCAMPO — segmentHA_T1.sh
# =============================================================================
echo ""
echo "--- PASO 2: segmentHA_T1.sh ---"
echo "Inicio: $(date)"

segmentHA_T1.sh "$ID_PACIENTE" "$SUBJECTS_DIR"

echo "Segmentación hipocampal finalizada: $(date)"
echo "Resultados en:"
echo "  ${SUBJECTS_DIR}/${ID_PACIENTE}/mri/lh.hippoSfVolumes-T1.v22.txt"
echo "  ${SUBJECTS_DIR}/${ID_PACIENTE}/mri/rh.hippoSfVolumes-T1.v22.txt"

# =============================================================================
# PASO 3: LESIONES DE SUSTANCIA BLANCA — SAMSEG con FLAIR
# =============================================================================
echo ""
echo "--- PASO 3: Lesiones de sustancia blanca (SAMSEG) ---"
echo "Inicio: $(date)"

cd "$FLAIR_DIR"

# 3a. Convertir DICOM FLAIR a NIfTI
echo "Convirtiendo FLAIR DICOM a NIfTI..."
~/bin/dcm2niix -z y -o . .

# Identificar el archivo .nii.gz generado.
# Si hay más de uno, el pipeline se detiene para evitar ambigüedad.
FLAIR_NIIGZ_COUNT=$(ls *.nii.gz 2>/dev/null | wc -l)
if [ "$FLAIR_NIIGZ_COUNT" -eq 0 ]; then
    echo "ERROR: dcm2niix no generó ningún archivo .nii.gz en $FLAIR_DIR. Abortando."
    exit 1
elif [ "$FLAIR_NIIGZ_COUNT" -gt 1 ]; then
    echo "ERROR: Se encontraron $FLAIR_NIIGZ_COUNT archivos .nii.gz en $FLAIR_DIR."
    echo "Archivos encontrados:"
    ls *.nii.gz
    echo "Editar el script para especificar manualmente la variable FLAIR_NII. Abortando."
    exit 1
fi

FLAIR_NII=$(ls *.nii.gz)
echo "Archivo FLAIR detectado: $FLAIR_NII"

# 3b. Convertir NIfTI a MGZ
mri_convert "$FLAIR_NII" flair.mgz
echo "Conversión a flair.mgz completada."

# 3c. Corregistro FLAIR → T1
echo "Corregistrando FLAIR a espacio T1..."
mri_coreg \
    --mov flair.mgz \
    --ref "${SUBJECTS_DIR}/${ID_PACIENTE}/mri/orig.mgz" \
    --reg flair2t1.lta

mri_vol2vol \
    --mov flair.mgz \
    --targ "${SUBJECTS_DIR}/${ID_PACIENTE}/mri/orig.mgz" \
    --reg flair2t1.lta \
    --o flair_in_t1.mgz \
    --interp nearest

echo "Corregistro completado."

# 3d. run_samseg — Volumen total de lesiones (genera samseg.stats)
echo "Ejecutando run_samseg para volumen total de lesiones..."
run_samseg \
    --input "${SUBJECTS_DIR}/${ID_PACIENTE}/mri/orig.mgz" flair_in_t1.mgz \
    --lesion \
    --output samseg_output

echo "Volumen total de lesiones guardado en: \
	${FLAIR_DIR}/samseg_output/samseg.stats"

# 3e. samseg — Mapa de lesiones con posteriors
echo "Ejecutando samseg para mapa de lesiones (posteriors)..."
samseg \
    --i flair_in_t1.mgz \
    --o samseg_output \
    --lesion \
    --save-posteriors

echo "Mapa de lesiones generado."

# 3f. Volumen de lesiones por región anatómica
echo "Calculando volumen de lesiones por región..."
mri_segstats \
    --seg "${SUBJECTS_DIR}/${ID_PACIENTE}/mri/aseg.mgz" \
    --mask samseg_output/posteriors/Lesions.mgz \
    --sum lesiones_por_region.txt

echo "Lesiones por región guardadas en: ${FLAIR_DIR}/lesiones_por_region.txt"

# =============================================================================
# RESUMEN FINAL
# =============================================================================
echo ""
echo "======================================================="
echo " Pipeline completado para: $ID_PACIENTE"
echo " Fecha/Hora: $(date)"
echo ""
echo " Archivos de resultados clave:"
echo "   [1] Volumetría global  : ${SUBJECTS_DIR}/${ID_PACIENTE}/stats/volumenes.txt"
echo "   [2] Hipocampo izq      : ${SUBJECTS_DIR}/${ID_PACIENTE}/mri/lh.hippoSfVolumes-T1.v22.txt"
echo "   [3] Hipocampo der      : ${SUBJECTS_DIR}/${ID_PACIENTE}/mri/rh.hippoSfVolumes-T1.v22.txt"
echo "   [4] Vol. total lesiones: ${FLAIR_DIR}/samseg_output/samseg.stats"
echo "   [5] Lesiones por región: ${FLAIR_DIR}/lesiones_por_region.txt"
echo "======================================================="
