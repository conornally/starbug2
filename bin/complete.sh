# file: starbug2

## need to make this local at some point
PARAMLIST=" PARAM\
TEST\ 
VERBOSE\ 
NULLVAL\ 
PSFDIR\ 
SIGSKY\ 
SIGSRC\ 
BOX_SIZE\ 
FILTER_SIZE\ 
SHARP_LO\ 
SHARP_HI\ 
ROUND_LO\ 
ROUND_HI\ 
APPHOT_R\ 
#APPHOT_R1\ 
#APPHOT_R2\ 
SKY_RIN\ 
SKY_ROUT\ 
ERROR_CUT\ 
SHARP_HI_SIG\ 
SHARP_LO_SIG\ 
ROUND_HI_SIG\ 
ROUND_LO_SIG\ 
BGD_R\ 
AP_FILE\ 
BGD_FILE\ 
CRIT_SEP\ 
MATCH_THRESH\ 
MATCH_COLS\ 
RM_MATCH\ 
NUMBER_ARTIFICIAL_STARS\ 
SUBIMAGE_SIZE\ 
MIN_FLUX\ 
MAX_FLUX\ 
SEPARATION_THRESH\ 
" 
_starbug2 ()
{
    COMPREPLY=()
    local cur=${COMP_WORDS[COMP_CWORD]}
    local prev=${COMP_WORDS[COMP_CWORD-1]}
    local xpat

    case "$prev" in 
      -p|--param)
          xpat="!*.param";;
      -b|--bgdfile) 
          xpat="!*-bgd.fits";;
      -d|--apfile) 
          xpat="!*-ap.fits";;

      ### THIS ONE STILL NEEDS SOME WORK
      -s|--set) 
          COMPREPLY=($(compgen -o nospace -W "$PARAMLIST" -- $cur))
          if [ ${#COMPREPLY[@]} == 1 ]
          then
              COMPREPLY=${COMPREPLY}"="
          fi
          return 0;;
      esac

    if [  $xpat ] 
    then
        COMPREPLY=($(compgen -o plusdirs -f -X "$xpat" -- $cur))
        return 0
    fi

    case "$cur" in 
    --*)
        COMPREPLY=( $( compgen -W ' --artific \
                                    --background \
                                    --bgdfile \
                                    --clean \
                                    --apfile \
                                    --detect \
                                    --find \
                                    --help \
                                    --output \
                                    --param \
                                    --photom \
                                    --set \
                                    --verbose \
                                    --generate-psf \
                                    --local-param \
                                    --generate-region \
                                    --clean-table \
                                    --version' -- $cur ) );;
    -*)
        COMPREPLY=($(compgen -W '-A -B -C -D -f -h -P -v \
                                -b -d -o -p -s ' -- $cur));;

    *) 
        xpat='!*.fits'
        COMPREPLY=($(compgen -o plusdirs -f -X "$xpat" -- $cur));;
    esac

    return 0
}

complete -o nospace -o filenames -F _starbug2 starbug2
