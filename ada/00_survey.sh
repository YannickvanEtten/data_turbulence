#!/bin/bash
# ===========================================================================
# ada/00_survey.sh -- WHAT DO I ALREADY HAVE?
#
# READ-ONLY. Creates nothing, changes nothing, deletes nothing. Run it before
# any setup so we know what is already on the cluster and do not clobber it.
#
# Run on the LOGIN node (it is light and finishes in seconds):
#     bash ada/00_survey.sh
#
# Or, if the repo is not on ADA yet, paste the whole file into a heredoc:
#     cat > survey.sh << 'ENDOFFILE'
#     ...paste...
#     ENDOFFILE
#     bash survey.sh
# ===========================================================================

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence

hr() { printf '\n%s\n' "-------------------------------------------------------------"; }
ok()   { printf '  [ OK ]   %s\n' "$1"; }
miss() { printf '  [MISSING] %s\n' "$1"; }
note() { printf '           %s\n' "$1"; }

echo "ADA survey -- $(date -Is)"
echo "user $(whoami) on $(hostname)"

hr; echo "1. THE PROJECT TREE"
if [ -d "$BASE" ]; then
  ok "$BASE exists"
  [ -w "$BASE" ] && ok "and is writable" || miss "but is NOT writable -- stop here, ask the group"
  echo
  echo "  what is in it now:"
  ls -la "$BASE" 2>/dev/null | tail -n +2 | sed 's/^/    /'
  echo
  echo "  space:"
  df -h "$BASE" 2>/dev/null | sed 's/^/    /'
else
  miss "$BASE does not exist"
  note "check the group share: ls -la /scistor/SBE-EDS-ClimateKoopman/"
fi

echo
echo "  subdirectories the pilot expects:"
for d in raw raw/north_atlantic raw/global derived calibration results logs; do
  if [ -d "$BASE/$d" ]; then
    n=$(find "$BASE/$d" -maxdepth 1 -type f 2>/dev/null | wc -l)
    ok "$d  ($n files at top level)"
  else
    miss "$d  (01_setup.sh creates it)"
  fi
done
# data_turbulence is deliberately NOT in that list: 01_setup.sh does not create
# the code, it requires it. Getting the repo onto ADA is a manual step (git
# clone, or MobaXterm SFTP) because it is a decision, not a default.

hr; echo "2. THE CODE"
if [ -d "$REPO/.git" ]; then
  ok "git repo at $REPO"
  note "branch:  $(git -C "$REPO" rev-parse --abbrev-ref HEAD 2>/dev/null)"
  note "commit:  $(git -C "$REPO" log -1 --format='%h %ad %s' --date=short 2>/dev/null)"
  note "dirty:   $(git -C "$REPO" status --porcelain 2>/dev/null | wc -l) modified file(s)"
elif [ -d "$REPO" ]; then
  ok "$REPO exists but is NOT a git checkout (copied by hand?)"
  note "updating it later will mean copying again rather than 'git pull'"
else
  miss "no code at $REPO -- see step 2 of the setup guide"
fi
echo "  key files:"
for f in pixi.toml 1_download_hpc.py download_plan.py 2_diagnostics.py 3_pipeline.py calibration.py; do
  [ -f "$REPO/$f" ] && ok "$f" || miss "$f"
done
[ -d "$REPO/tests" ] && ok "tests/ ($(ls "$REPO"/tests/*.py 2>/dev/null | wc -l) files)" || miss "tests/"
[ -d "$REPO/jobs" ] && ok "jobs/ ($(ls "$REPO"/jobs/*.sbatch 2>/dev/null | wc -l) sbatch)" || miss "jobs/"

hr; echo "3. THE ENVIRONMENT"
if module load 2025 >/dev/null 2>&1; then ok "module load 2025"; else miss "module load 2025 FAILED"; fi
if module load pixi >/dev/null 2>&1; then
  ok "module load pixi  ->  $(pixi --version 2>&1)"
else
  miss "module load pixi FAILED (note: lowercase 'pixi', not 'Pixi')"
fi

echo
echo "  existing pixi environments:"
found_env=0
for d in "$REPO" "$HOME/env-test" "$BASE/env"; do
  if [ -f "$d/pixi.toml" ]; then
    found_env=1
    ok "pixi.toml at $d"
    if [ -d "$d/.pixi" ]; then
      note "BUILT ($(du -sh "$d/.pixi" 2>/dev/null | cut -f1))"
      rv=$(cd "$d" && pixi run python -c "import importlib.metadata as m;print(m.version('rojak-cat'))" 2>/dev/null)
      [ -n "$rv" ] && note "rojak-cat $rv" || note "rojak NOT importable in it"
    else
      note "not built yet (no .pixi/)"
    fi
    note "rojak rev pinned: $(grep -o 'rev = \"[^\"]*\"' "$d/pixi.toml" 2>/dev/null | head -1)"
  fi
done
[ "$found_env" = 0 ] && miss "no pixi.toml found anywhere -- 01_setup.sh will build one"
echo
note "NOTE: ~/env-test was a THROWAWAY validation env from July. The pilot"
note "      builds a fresh one inside the repo. Leave env-test alone; it costs"
note "      nothing and is a working fallback if the new build misbehaves."

hr; echo "4. CDS CREDENTIALS"
if [ -f ~/.cdsapirc ]; then
  ok "~/.cdsapirc present, mode $(stat -c '%a' ~/.cdsapirc)"
  [ "$(stat -c '%a' ~/.cdsapirc)" = "600" ] || note "mode should be 600:  chmod 600 ~/.cdsapirc"
  u=$(grep -i '^url' ~/.cdsapirc | head -1)
  note "$u"
  case "$u" in
    *api/v2*) note "!! OLD CDS format. The new API needs 'url: https://cds.climate.copernicus.eu/api'" ;;
    *)        note "looks like the new format" ;;
  esac
  note "(the key itself is deliberately not printed)"
else
  miss "~/.cdsapirc -- every download will fail without it"
fi

hr; echo "5. SCHEDULER"
sinfo -s 2>/dev/null | head -8 | sed 's/^/    /'
echo "  your queued/running jobs:"
squeue --me 2>/dev/null | sed 's/^/    /'
echo "  recent job history:"
sacct -u "$USER" -S "$(date -d '30 days ago' +%F 2>/dev/null || echo 2026-07-01)" \
      --format=JobID%12,JobName%18,State%12,Elapsed%10 2>/dev/null | head -12 | sed 's/^/    /'

hr; echo "SURVEY DONE"
echo "Paste this output back and we will pick up from exactly where you are."
