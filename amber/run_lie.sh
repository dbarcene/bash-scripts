#!/usr/bin/env bash
#
#SBATCH --job-name=lie
#SBATCH --partition=gpu
#SBATCH --nodelist=gpu01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

scl enable devtoolset-9 bash
set -euo pipefail

# ==========================================================
#   Single-trajectory LIE (stLIE) + RMSD/RMSF + contacts
#   ligand–residues (5 A) + residue-level summary via PDB
# ==========================================================

# -------- LIE parameters --------
alpha=0.18
beta=0.50
gamma=0.0

use_water_terms=0
alpha_w=0.00
beta_w=0.00

# -------- MD inputs --------
PRMTOP_COMPLEX="complex_solv.top"
TRAJ_COMPLEX="pol2_1_dyn.nc"
TRAJ_COMPLEX1="pol2_2_dyn.nc"

# -------- PDB for residue name mapping --------
# Debe ser consistente con la topologia (mismos numeros de residuo)
PDB_FILE="limpio.pdb"

# -------- Ligand residue ID (VERY IMPORTANT) --------
# Ajusta este numero al residuo del ligando en tu sistema (segun tleap/top)
LIG_RESID=353

# -------- Masks --------
LIG_MASK=":UNL"
REC_MASK="(:1-353)&!(:WAT,Na+,Cl-,UNL)"        # protein only, no water/ions/ligand
WAT_MASK=":WAT,Na+,Cl-"
TOT_MASK="(:1-353,UNL)&!(:WAT,Na+,Cl-)"       # protein + ligand, no solvent

# -------- Trajectory range --------
start=1
end=10000
stride=10

# Integration time step (ps)
dt_ps=0.004

OUTDIR="stLIE_out"
mkdir -p "$OUTDIR"

# ==========================================================
#   1) CPPTRAJ: LIE, RMSD, RMSF, ligand–residue contacts
# ==========================================================
cat > "$OUTDIR/analysis.in" <<EOF
parm $PRMTOP_COMPLEX
trajin $TRAJ_COMPLEX  $start $end $stride
trajin $TRAJ_COMPLEX1 $start $end $stride

# ---------- LIE ----------
lie $LIG_MASK $REC_MASK out $OUTDIR/LIG_REC.dat diel 80
$( [ "$use_water_terms" -eq 1 ] && echo "lie $LIG_MASK $WAT_MASK out $OUTDIR/LIG_WAT.dat diel 80" )

# ---------- RMSD ----------
rms RMS_tot   first $TOT_MASK  out $OUTDIR/RMSD_total.dat
rms RMS_lig   first $LIG_MASK  out $OUTDIR/RMSD_UNL.dat
rms RMS_prot  first $REC_MASK  out $OUTDIR/RMSD_protein.dat

# ---------- RMSF ----------
rmsf RMSF_prot $REC_MASK byres out $OUTDIR/RMSF_protein.dat
rmsf RMSF_lig  $LIG_MASK       out $OUTDIR/RMSF_UNL.dat

# ---------- Ligand–residue contacts (<= 5 A) ----------
nativecontacts name LIG_RES \
  :UNL $REC_MASK \
  distance 5.0 byresidue \
  out $OUTDIR/LIG_contacts_5A.dat \
  resout $OUTDIR/LIG_contacts_residue.dat

run
EOF

echo ">>> Running cpptraj..."
cpptraj -i "$OUTDIR/analysis.in" > "$OUTDIR/cpptraj.log" 2>&1

# ==========================================================
#   2) LIE, RMSD, RMSF statistics (AWK)
# ==========================================================
total_frames=$(grep -v '^#' "$OUTDIR/RMSD_total.dat" | awk 'NF>1' | wc -l)

time_ns=$(awk -v n="$total_frames" -v dt="$dt_ps" -v stride="$stride" '
  BEGIN{
    if (n>1) printf("%.3f", (n-1)*dt*stride/1000.0);
    else     printf("0.000");
  }')

# --- AWK: LIE ---
calc_mean_sd_LIE() {
  local file="$1"
  awk '
    BEGIN{vdw_col=0; ele_col=0}
    $1=="#Frame"||$1=="Frame"{
      for(i=1;i<=NF;i++){
        if(toupper($i)~/VDW/)  vdw_col=i
        if(toupper($i)~/ELEC/) ele_col=i
      }
      next
    }
    $0!~/^#/{
      if(vdw_col==0||ele_col==0){vdw_col=2;ele_col=3}
      vdw=$vdw_col+0; ele=$ele_col+0
      vdw_sum+=vdw; ele_sum+=ele
      vdw2_sum+=vdw*vdw; ele2_sum+=ele*ele
      n++
    }
    END{
      if(n>1){
        vdw_mean=vdw_sum/n; ele_mean=ele_sum/n
        vdw_sd=sqrt((vdw2_sum - n*vdw_mean^2)/(n-1))
        ele_sd=sqrt((ele2_sum - n*ele_mean^2)/(n-1))
        printf("%.6f %.6f %.6f %.6f\n", vdw_mean, ele_mean, vdw_sd, ele_sd)
      } else print "nan nan nan nan"
    }' "$file"
}

# --- AWK: RMSD ---
calc_mean_sd_RMSD() {
  local file="$1"
  awk '
    $1=="#Frame"||$1=="Frame"{next}
    $0!~/^#/{
      if(NF<2) next
      x=$2+0
      sum+=x
      sum2+=x*x
      n++
    }
    END{
      if(n>1){
        mean=sum/n
        sd=sqrt((sum2 - n*mean^2)/(n-1))
        printf("%.6f %.6f\n", mean, sd)
      } else print "nan nan"
    }' "$file"
}

# --- AWK: RMSF ---
calc_mean_sd_RMSF() {
  local file="$1"
  awk '
    $1=="#Res"||$1=="#Atom"||$1=="Res"||$1=="Atom"{next}
    $0!~/^#/{
      if(NF<2) next
      x=$2+0
      sum+=x
      sum2+=x*x
      n++
    }
    END{
      if(n>1){
        mean=sum/n
        sd=sqrt((sum2 - n*mean^2)/(n-1))
        printf("%.6f %.6f\n", mean, sd)
      } else print "nan nan"
    }' "$file"
}

# --- LIE: averages ---
read E_vdw_LR E_ele_LR SD_vdw_LR SD_ele_LR < <(calc_mean_sd_LIE "$OUTDIR/LIG_REC.dat")

if [ "$use_water_terms" -eq 1 ]; then
  read E_vdw_LW E_ele_LW SD_vdw_LW SD_ele_LW < <(calc_mean_sd_LIE "$OUTDIR/LIG_WAT.dat")
else
  E_vdw_LW=0.0;  E_ele_LW=0.0
  SD_vdw_LW=0.0; SD_ele_LW=0.0
fi

# --- Delta G bind and its SD ---
read DG_bind SD_DG < <(awk -v a="$alpha" -v b="$beta" -v g="$gamma" \
                           -v aw="$alpha_w" -v bw="$beta_w" \
                           -v vdwLR="$E_vdw_LR" -v eleLR="$E_ele_LR" \
                           -v vdwLW="$E_vdw_LW" -v eleLW="$E_ele_LW" \
                           -v sd_vdw="$SD_vdw_LR" -v sd_ele="$SD_ele_LR" '
  function s(x){return (x==""||x=="nan"||x=="NaN")?0.0:x+0.0}
  BEGIN{
    DG    = a*s(vdwLR) + b*s(eleLR) + aw*s(vdwLW) + bw*s(eleLW) + g
    SD_DG = sqrt( (a*s(sd_vdw))^2 + (b*s(sd_ele))^2 )
    printf("%.6f %.6f\n", DG, SD_DG)
  }')

# --- RMSD stats ---
read RMS_tot_mean   RMS_tot_sd   < <(calc_mean_sd_RMSD "$OUTDIR/RMSD_total.dat")
read RMS_lig_mean   RMS_lig_sd   < <(calc_mean_sd_RMSD "$OUTDIR/RMSD_UNL.dat")
read RMS_prot_mean  RMS_prot_sd  < <(calc_mean_sd_RMSD "$OUTDIR/RMSD_protein.dat")

# --- RMSF stats ---
read RMSF_prot_mean RMSF_prot_sd < <(calc_mean_sd_RMSF "$OUTDIR/RMSF_protein.dat")
read RMSF_lig_mean  RMSF_lig_sd  < <(calc_mean_sd_RMSF "$OUTDIR/RMSF_UNL.dat")

# ==========================================================
#   3) Residue-level contact analysis using PDB (Python)
# ==========================================================
CONTACT_SUMMARY="$OUTDIR/LIG_contacts_summary.dat"

export PDB_FILE
export OUTDIR
export CONTACT_SUMMARY
export LIG_RESID

python3 << 'PY'
import os

pdb_file = os.environ.get("PDB_FILE", "")
outdir = os.environ.get("OUTDIR", ".")
out_summary = os.environ.get("CONTACT_SUMMARY", os.path.join(outdir, "LIG_contacts_summary.dat"))
res_contacts_file = os.path.join(outdir, "LIG_contacts_residue.dat")

try:
    lig_resid = int(os.environ.get("LIG_RESID", "0"))
except ValueError:
    lig_resid = 0

res_info = {}      # resid -> (resname, chain)

# --- 1) Read PDB to map residue IDs to names and chains ---
if pdb_file and os.path.isfile(pdb_file):
    with open(pdb_file, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            resname = line[17:20].strip()
            chain   = line[21].strip()
            try:
                resid = int(line[22:26])
            except ValueError:
                continue
            if resid not in res_info:
                res_info[resid] = (resname, chain)
else:
    print("WARNING: PDB file not found or not set:", pdb_file)

# --- 2) Read residue-level contacts from nativecontacts (resout) ---
contact_acc = {}  # prot_resid -> total_fraction

if os.path.isfile(res_contacts_file):
    with open(res_contacts_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                r1 = int(parts[0])
                r2 = int(parts[1])
                frac = float(parts[2])
            except ValueError:
                continue

            # Decide which residue is ligand and which is protein using LIG_RESID
            if r1 == lig_resid and r2 != lig_resid:
                prot = r2
            elif r2 == lig_resid and r1 != lig_resid:
                prot = r1
            else:
                continue

            contact_acc[prot] = contact_acc.get(prot, 0.0) + frac
else:
    print("WARNING: contact file not found:", res_contacts_file)

# --- 3) Sort and write summary ---
rows = []
for resid, frac in contact_acc.items():
    resname, chain = res_info.get(resid, ("UNK", ""))
    rows.append((frac, resid, resname, chain))

rows.sort(key=lambda x: x[0], reverse=True)

with open(out_summary, "w") as out:
    out.write("#Resid  ResName Chain  ContactFraction_le_5A_with_LIG\n")
    for frac, resid, resname, chain in rows:
        out.write(f"{resid:5d}  {resname:>6s}  {chain:1s}  {frac:.4f}\n")

print("Wrote contact summary to", out_summary)
PY

# ==========================================================
#   4) Final report
# ==========================================================
{
  echo "==============================================="
  echo " Single-trajectory LIE (stLIE) — Results"
  echo "==============================================="
  echo "Trajectories: $TRAJ_COMPLEX + $TRAJ_COMPLEX1"
  echo "Masks:"
  echo "  LIG  = $LIG_MASK  (resid $LIG_RESID)"
  echo "  REC  = $REC_MASK"
  echo "  TOT  = $TOT_MASK"
  echo "  WAT  = $WAT_MASK"
  echo "Frames analyzed (RMSD_total): $total_frames"
  echo "Effective time: $time_ns ns  (dt = $dt_ps ps, stride = $stride)"
  echo
  echo "LIE (kcal/mol):"
  echo "  <VvdW(LIG-REC)> = $E_vdw_LR ± $SD_vdw_LR"
  echo "  <Vele(LIG-REC)> = $E_ele_LR ± $SD_ele_LR"
  if [ "$use_water_terms" -eq 1 ]; then
    echo "  <VvdW(LIG-WAT)> = $E_vdw_LW ± $SD_vdw_LW"
    echo "  <Vele(LIG-WAT)> = $E_ele_LW ± $SD_ele_LW"
  fi
  echo
  echo "LIE parameters:"
  echo "  alpha = $alpha, beta = $beta, gamma = $gamma"
  if [ "$use_water_terms" -eq 1 ]; then
    echo "  alpha_w = $alpha_w, beta_w = $beta_w"
  fi
  echo "-----------------------------------------------"
  echo "DG_bind (kcal/mol) = $DG_bind ± $SD_DG"
  echo "-----------------------------------------------"
  echo
  echo "RMSD (A):"
  echo "  RMSD total   (protein + ligand, no solvent) = $RMS_tot_mean  ± $RMS_tot_sd"
  echo "  RMSD ligand                              = $RMS_lig_mean  ± $RMS_lig_sd"
  echo "  RMSD protein (REC_MASK)                  = $RMS_prot_mean ± $RMS_prot_sd"
  echo
  echo "RMSF (A):"
  echo "  RMSF protein (by residue)                = $RMSF_prot_mean ± $RMSF_prot_sd"
  echo "    -> Full profile: $OUTDIR/RMSF_protein.dat"
  echo "  RMSF ligand  (by atom)                   = $RMSF_lig_mean  ± $RMSF_lig_sd"
  echo "    -> Full profile: $OUTDIR/RMSF_UNL.dat"
  echo
  echo "Ligand–residue contacts (<= 5 A):"
  echo "  cpptraj resout file: $OUTDIR/LIG_contacts_residue.dat"
  echo "  Residue-level summary (using PDB: $PDB_FILE):"
  echo "    -> $CONTACT_SUMMARY"
  echo "==============================================="
} | tee "$OUTDIR/stLIE_report.txt"
