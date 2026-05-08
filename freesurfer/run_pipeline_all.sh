#!/bin/bash
# File              : run_pipeline_all.sh
# Author            : David Barcene <dbarcene@indicasat.org>
# Date              : 29.04.2026
# Last Modified Date: 30.04.2026
# Last Modified By  : David Barcene <dbarcene@indicasat.org.pa>
#SBATCH --job-name=pipeline_all
#SBATCH --output=/data/aberraondo/FS_LOGS/fs_%A_%a.out
#SBATCH --error=/data/aberraondo/FS_LOGS/fs_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --array=2-117%4
#
# USO:
#   mkdir -p /home/aberraondo/FS_LOGS       # EJECUTAR ANTES de sbatch
#   sbatch run_pipeline_all.sh


# =============================================================================
# 1. Configuración del Entorno
# =============================================================================
module purge
module load freesurfer/8.1.0-1
# El modulefile ya ejecuta: export FREESURFER_HOME=...; source SetUpFreeSurfer.sh
# El modulefile ya setea:   FS_LICENSE, SUBJECTS_DIR, FUNCTIONALS_DIR
# Es necesario redefinir SUBJECTS_DIR si se quiere una ruta distinta.

# Cambiar SUBJECTS_DIR a /data/aberraondo/RESULTADOS u otra carpeta
export SUBJECTS_DIR="/data/aberraondo/RESULTADOS"
mkdir -p "$SUBJECTS_DIR"

# =============================================================================
# 2. Mapeo de Pacientes
# =============================================================================
PATIENT_DIRS=(/data/aberraondo/VISITA_1/03-*)

# Validar que el task ID está dentro del rango del array real
if [ "$SLURM_ARRAY_TASK_ID" -ge "${#PATIENT_DIRS[@]}" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID \
	    excede el número de pacientes (${#PATIENT_DIRS[@]}). Abortando."
    exit 1
fi

CURRENT_PATH="${PATIENT_DIRS[$SLURM_ARRAY_TASK_ID]}"
PACIENTE=$(basename "$CURRENT_PATH")

DICOM_BASE="$CURRENT_PATH/DICOM/00000001"
T1_DIR="$DICOM_BASE/00000002"
FLAIR_DIR="$DICOM_BASE/00000003"

echo "======================================================="
echo " Procesando Paciente : $PACIENTE"
echo " SLURM Task ID       : $SLURM_ARRAY_TASK_ID"
echo " Ruta T1             : $T1_DIR"
echo " Ruta FLAIR          : $FLAIR_DIR"
echo " SUBJECTS_DIR        : $SUBJECTS_DIR"
echo " Fecha/Hora          : $(date)"
echo "======================================================="

# =============================================================================
# 3. Volumetría del Hipocampo (Recon-all)
# =============================================================================
echo "--- PASO 1: recon-all --- $(date)"

recon-all -all \
    -i "$T1_DIR/00000001" \
    -s "$PACIENTE" \
    -sd "$SUBJECTS_DIR" \
    -threads 20

# Generar tabla de volúmenes globales
cd "$SUBJECTS_DIR/$PACIENTE/stats"
asegstats2table --subjects "$PACIENTE" --meas volume --tablefile volumenes.txt

# =============================================================================
# 4. Segmentación del Hipocampo
# =============================================================================
echo "--- PASO 2: segmentHA_T1.sh --- $(date)"

cd "$HOME"
segmentHA_T1.sh "$PACIENTE" "$SUBJECTS_DIR"

# =============================================================================
# 5. Lesiones de Sustancia Blanca (SAMSEG)
# =============================================================================
echo "--- PASO 3: SAMSEG --- $(date)"

cd "$FLAIR_DIR"
~/bin/dcm2niix -z y -o . .

# FIX #4: Reemplazar wildcard *.nii.gz por detección explícita.
# En un batch sin supervisión un wildcard ambiguo produce fallos silenciosos.
FLAIR_NIIGZ_COUNT=$(ls *.nii.gz 2>/dev/null | wc -l)
if [ "$FLAIR_NIIGZ_COUNT" -eq 0 ]; then
    echo "ERROR [$PACIENTE]: dcm2niix no generó ningún .nii.gz en $FLAIR_DIR. Abortando."
    exit 1
elif [ "$FLAIR_NIIGZ_COUNT" -gt 1 ]; then
    echo "ERROR [$PACIENTE]: $FLAIR_NIIGZ_COUNT archivos .nii.gz en $FLAIR_DIR — ambigüedad. Abortando."
    ls *.nii.gz
    exit 1
fi
FLAIR_NII=$(ls *.nii.gz)
echo "Archivo FLAIR: $FLAIR_NII"

mri_convert "$FLAIR_NII" flair.mgz

# Corregistro FLAIR → T1
mri_coreg \
    --mov flair.mgz \
    --ref "$SUBJECTS_DIR/$PACIENTE/mri/orig.mgz" \
    --reg flair2t1.lta

mri_vol2vol \
    --mov flair.mgz \
    --targ "$SUBJECTS_DIR/$PACIENTE/mri/orig.mgz" \
    --reg flair2t1.lta \
    --o flair_in_t1.mgz \
    --interp nearest

# FIX #3: Añadir run_samseg para generar samseg.stats (volumen total de lesiones).
# Estaba ausente en el script original; sin este paso samseg.stats no se crea.
echo "Ejecutando run_samseg (volumen total)..."
run_samseg \
    --input "$SUBJECTS_DIR/$PACIENTE/mri/orig.mgz" flair_in_t1.mgz \
    --lesion \
    --output samseg_output

# Mapa de lesiones con posteriors
echo "Ejecutando samseg (mapa de posteriors)..."
samseg \
    --i flair_in_t1.mgz \
    --o samseg_output \
    --lesion \
    --save-posteriors

# Volumen de lesiones por región anatómica
mri_segstats \
    --seg "$SUBJECTS_DIR/$PACIENTE/mri/aseg.mgz" \
    --mask samseg_output/posteriors/Lesions.mgz \
    --sum lesiones_por_region.txt

echo "======================================================="
echo " Pipeline completado: $PACIENTE — $(date)"
echo "======================================================="
