#!/bin/bash
# Set up the villa spiral fitter under WSL Ubuntu (Linux build of the native vc_spiral extension, torch cu128).
W="$(cd "$(dirname "$0")" && pwd)"   # this scripts/ folder
set -e
export PATH="$HOME/.local/bin:$PATH"
command -v uv >/dev/null || curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"
mkdir -p ${R3:-$HOME/r3}
rsync -a --delete --exclude .venv "${VILLA:?set VILLA to your villa checkout}"/spiral-fitting/ ${R3:-$HOME/r3}/spiral-fitting/
cd ${R3:-$HOME/r3}/spiral-fitting
uv sync 2>&1 | tail -20
.venv/bin/python -c "import vc_spiral.spiral_sampling, vc_spiral.track_store, torch, triton; print('NATIVE OK', torch.__version__, torch.cuda.is_available(), triton.__version__)"
echo WSL_SETUP_DONE
