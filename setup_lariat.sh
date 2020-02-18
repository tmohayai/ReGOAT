# To setup a new LArIATSoft development area, follow the instructions on https://cdcvs.fnal.gov/redmine/projects/goat/wiki/Setting_up_Offline_Software
# note: with a recent upgrade to dunegpvms, only sl7 would work for setting up LArIATSoft. 
# Run this set up script after each login once a LArIATSoft development area is set up

cd to_your_lariat_directory/
source /grid/fermiapp/larsoft/products/setup
export PRODUCTS=/grid/fermiapp/products/lariat:${PRODUCTS}:/grid/fermiapp/products/common/db
setup git
setup gitflow
setup mrb
setup ninja v1_6_0b

source localProducts_lariatsoft_develop_e14_prof/setup
mrbsetenv
setup lariatsoft v06_62_00 -q e14:prof
setup gallery v1_05_03 -q e14:nu:prof
