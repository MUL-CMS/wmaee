#!/bin/bash
# >>> mamba initialize >>>
# !! Contents within this block are managed by 'mamba shell init' !!
export MAMBA_EXE='/home/dholec/software/miniconda3/bin/mamba';
export MAMBA_ROOT_PREFIX='/home/dholec/software/miniconda3';
__mamba_setup="$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__mamba_setup"
else
    alias mamba="$MAMBA_EXE"  # Fallback on help from mamba activate
fi
unset __mamba_setup
# <<< mamba initialize <<<

mamba activate grace

export LD_PRELOAD=$CONDA_PREFIX/lib/libstdc++.so.6
python calc_script.py
