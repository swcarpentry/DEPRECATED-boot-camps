# /etc/skel/.bashrc
#
# This file is sourced by all *interactive* bash shells on startup,
# including some apparently interactive shells such as scp and rcp
# that can't tolerate any output.  So make sure this doesn't display
# anything or bad things will happen !

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
# Test for an interactive shell.  There is no need to set anything
# past this point for scp and rcp, and it's important to refrain from
# outputting anything in those cases.
if [[ $- != *i* ]] ; then
	# Shell is non-interactive.  Be done now!
	return
fi


# Put your fun stuff here!
export EDITOR="nano"

#
# Version control current branch function
#
function get_curr_branch {
  local dir="$PWD"
  local vcs
  local nick
  while [[ "$dir" != "/" ]]; do
    for vcs in git hg svn bzr; do
      if [[ -d "$dir/.$vcs" ]] && hash "$vcs" &>/dev/null; then
        case "$vcs" in
          git) __git_ps1 "${1:- %s}"; return;;
          hg) nick=$(hg branch 2>/dev/null);;
          svn) nick=$(svn info 2>/dev/null\
                | grep -e '^Repository Root:'\
                | sed -e 's#.*/##');;
          bzr)
            local conf="${dir}/.bzr/branch/branch.conf" # normal branch
            [[ -f "$conf" ]] && nick=$(grep -E '^nickname =' "$conf" | cut -d' ' -f 3)
            conf="${dir}/.bzr/branch/location" # colo/lightweight branch
            [[ -z "$nick" ]] && [[ -f "$conf" ]] && nick="$(basename "$(< $conf)")"
            [[ -z "$nick" ]] && nick="$(basename "$(readlink -f "$dir")")";;
        esac
        [[ -n "$nick" ]] && printf "${1:- %s}" "$nick"
        return 0
      fi
    done
    dir="$(dirname "$dir")"
  done
}


# Git Friendly PS1
#PS1='\[\033[01;32m\]\u@\h\[\033[01;34m\] \w\[\033[01;31m\]`git branch 2>/dev/null | grep ^* | tr -d \*` \[\033[01;34m\]\$\[\033[00m\] '

# Version Control Friendly PS1
PS1='\[\033[01;32m\]\u@\h\[\033[01;34m\] \w\[\033[01;31m\]$(get_curr_branch '$2') \[\033[01;34m\]\$\[\033[00m\] '

# Hacky fix for less
export PAGER='most'
alias less='most'


# Enable custom input control
HISTSIZE=8128

#
# Alias definitions.
#
# Color aliases
alias grep='grep --color=auto'
alias fgrep='fgrep --color=auto'
alias egrep='egrep --color=auto'

# Correct silly ubuntu mistakes
alias ls='ls --color=auto -v'

# Some awesome aliases
alias scp-resume="rsync --partial -h --progress --rsh=ssh"
alias git-commit-count="git log --pretty=format:'' | wc -l"
alias pep8="pep8 --ignore=E501 -r"

# Enable programmable completion features (you don't need to enable
# this, if it's already enabled in /etc/bash.bashrc and /etc/profile
# sources /etc/bash.bashrc).
if [ -f /etc/bash_completion ] && ! shopt -oq posix; then
    . /etc/bash_completion
fi

