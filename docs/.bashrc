# shell config for use in Docker for Jekyll

export PS1="\e[0;32m\w > \e[0m"

export LS_OPTIONS='--color=auto'
export LS_COLORS='rs=0:di=01;34:ln=01;36:mh=00:pi=40;33:so=01;35:do=01;35:bd=40;33;01:cd=40;33;01:or=40;31;01:mi=00:su=37;41:sg=30;43:ca=30;41:tw=30;42:ow=34;42:st=37;44:ex=01;32:*.vcf=00;35:*.preprocessed.vcf=00;33:*.match.json=00;33:*.phenotype.json=00;33:*.report.html=00;33:'
alias ls='ls $LS_OPTIONS'
alias ll='ls $LS_OPTIONS -l'
alias l='ls $LS_OPTIONS -lA'

alias rm='rm -i'
alias cp='cp -i'
alias mv='mv -i'
