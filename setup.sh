source /grid/fermiapp/products/uboone/setup_uboone.sh; 
setup git v2_15_1; 
setup gitflow v1_11_0; 
setup mrb v1_16_02;
setup root v6_08_06g -q e14:nu:prof; 
setup cmake v2_8_8;


me="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $me

export SBNFITDIR=$me
export SBNFIT_LIBDIR=$me/build

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SBNFIT_LIBDIR




