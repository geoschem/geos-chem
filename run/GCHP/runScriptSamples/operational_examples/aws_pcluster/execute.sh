#!/bin/bash

if [[ $( cat /proc/sys/kernel/yama/ptrace_scope ) == "0" ]]; then
        echo "PTrace Correct for EFA"
else
        echo "PTrace Override"
        sysctl -w kernel.yama.ptrace_scope=0
fi

./gchp