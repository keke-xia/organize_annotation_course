#!/bin/bash
# Local R environment setup for VSCode Remote SSH (fixed)
set -euo pipefail

# ----- Config -----
WORKDIR="/data/users/kxia/organize_annotation_course"
LOGDIR="$WORKDIR/logs"
ENV_NAME="r-env"
R_USER_LIB="$WORKDIR/.Rlibs_renv"

mkdir -p "$LOGDIR" "$R_USER_LIB"

echo "[I] Local R env setup start"

# 0) Kill module R influence; only use conda's R
module purge || true

# 1) Load Anaconda and activate env
module load Anaconda3
source "$(conda info --base)/etc/profile.d/conda.sh"

if ! conda env list | grep -q "^${ENV_NAME}\b"; then
  echo "[I] Creating env '${ENV_NAME}' with r-base..."
  conda create -n "$ENV_NAME" -y conda-forge::r-base
fi
conda activate "$ENV_NAME"
echo "[✓] Activated env: $ENV_NAME"

# 2) Force R to use writable user library
export R_LIBS_USER="$R_USER_LIB"

# 3) Install packages into R_LIBS_USER
Rscript -e '
userlib <- Sys.getenv("R_LIBS_USER")
dir.create(userlib, showWarnings=FALSE, recursive=TRUE)
.libPaths(c(userlib, .libPaths()))  # put user lib first
repos <- c(CRAN="https://cloud.r-project.org")
pkgs <- c("languageserver","tidyverse","data.table","circlize","cowplot","reshape2")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need, repos=repos, lib=userlib, Ncpus=2)
cat("[I] .libPaths():\n"); print(.libPaths())
cat("[I] Installed versions:\n")
print(installed.packages()[pkgs, c("Package","Version"), drop=FALSE])
'

echo "[✓] R env ready at $ENV_NAME ; user library: $R_USER_LIB"
echo "[i] Next time use:"
echo "    screen -R R_terminal"
echo "    module purge && module load Anaconda3 && conda activate r-env"
echo "    R"

# ---------- Step 5: Usage instructions ----------
cat <<'EOT'

==============================================================
[✓] R environment setup complete!
Next steps to use R in VSCode:
--------------------------------------------------------------
1. On your local computer:
   - Open VSCode → Extensions
   - Install **R Editor Support** (`REditorSupport.r`)

2. In VSCode → Settings:
   - Search: `@ext:REditorSupport.r R: Always use Active Terminal`
   - ✅ Check the box to enable it.

3. Connect to HPC via Remote SSH in VSCode:
   - Open folder: /data/users/kxia/organize_annotation_course
   - Open an .R file (e.g. plot_div.R)

4. In the VSCode terminal:
   - Run:
        screen -R R_terminal
        conda activate r-env
        R
   - Then in the R console, run:
        install.packages("languageserver")

5. Now you can run R code directly in VSCode with **Ctrl + Enter**.
==============================================================

EOT