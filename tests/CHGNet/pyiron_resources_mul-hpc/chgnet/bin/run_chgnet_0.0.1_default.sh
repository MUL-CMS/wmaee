#!/bin/bash


# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/mulfs/home/p0817489/software/mambaforge/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/mulfs/home/p0817489/software/mambaforge/etc/profile.d/conda.sh" ]; then
        . "/mulfs/home/p0817489/software/mambaforge/etc/profile.d/conda.sh"
    else
        export PATH="/mulfs/home/p0817489/software/mambaforge/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "/mulfs/home/p0817489/software/mambaforge/etc/profile.d/mamba.sh" ]; then
    . "/mulfs/home/p0817489/software/mambaforge/etc/profile.d/mamba.sh"
fi
# <<< conda initialize <<<
conda activate chgnet

python calc_script.py

