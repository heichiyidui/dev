# Lines configured by zsh-newuser-install
HISTFILE=~/.histfile
HISTSIZE=1000
SAVEHIST=1000
setopt autocd beep
unsetopt extendedglob nomatch
bindkey -v
# End of lines configured by zsh-newuser-install

# The following lines were added by compinstall

zstyle ':completion:*' completer _complete _ignored
zstyle ':completion:*' max-errors 1
zstyle :compinstall filename '/home/kuang/.zshrc'

autoload -Uz compinit
compinit
# End of lines added by compinstall

autoload -U promptinit
promptinit
autoload -U colors && colors
PROMPT="%{$fg[red]%T%} %U%{$fg[green]%~%}%u > %{$reset_color%}"

export PATH=.:~/bin:$PATH
export LANG=en_GB

test -s ~/.alias && . ~/.alias || true

# an2WH2jW9Hw5

