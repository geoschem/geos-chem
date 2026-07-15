#!/bin/bash
41;330;0c
alias compile='source gchp.env; cd build; make -j; cp bin/gchp ..; cd ..'
alias run='./cleanRunDir.sh; sbatch gchp.run.mapl3'

alias log='emacs run.log -nw'
alias logc='cat run.log'
alias logh='head run.log -n 30'
alias logt='tail run.log -n 30'

alias all='emacs allPEs.log -nw'
alias allc='cat allPEs.log'
alias allh='head allPEs.log -n 30'
alias allt='tail allPEs.log -n 30'

alias pet='emacs PET0.ESMF_LogFile -nw'
alias petc='cat PET0.ESMF_LogFile'
alias peth='head PET0.ESMF_LogFile -n 30'
alias pett='tail PET0.ESMF_LogFile -n 30'

