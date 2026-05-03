#!/bin/bash
body() {
    IFS= read -r header
    printf '%s\n' "$header"
    "$@"
}
xjoin() {
    local f
    local srt="sort -k 1b,1"

    if [ "$#" -lt 2 ]; then
            echo "xjoin: need at least 2 files" >&2
            return 1
    elif [ "$#" -lt 3 ]; then
            join -t "	" <(tail -n+2 "$1" | $srt) <(tail -n+2 "$2" | $srt)
    else
            f=$1
            shift
            join -t "	" <(tail -n+2 "$f" | $srt) <(xjoin "$@")
    fi
}

xjoin_head() {
	local f
	if [ "$#" -lt 2 ]; then
		head -1 "$1" | cut -d '	' -f2-	
	elif [ "$#" -lt 3 ]; then
		head -1 "$1" | cut -d '	' -f2- | paste - <(head -1 "$2" | cut -d '	' -f2-)
	else
		f=$1
		shift
		head -1 "$f" | cut -d '	' -f2 | paste - <(xjoin_head "$@") 
	fi
}
cat <(head -1 "$1" | paste - <(xjoin_head ${@:2})) <(xjoin $* | sort -k 1n,1)
