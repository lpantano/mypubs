#nice options for bash scripts
set -o pipefail # trace ERR through pipes
set -o errtrace # trace ERR through 'time command' and other functions
set -o nounset ## set -u : exit the script if you try to use an uninitialised variable
set -o errexit ## set -e : exit the script if any statement returns a non-true return value
set -v
set -x
export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

razers3 -dr 0 -i 80 -rr 90 -f -o hsa_let_7a.sam hairpin.hsa.fa hsa_let_7a.fa
