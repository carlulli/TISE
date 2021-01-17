################################################################################
# This script can be used to execute the simulation, namely running the main.o
# like this:
# ./main.o [int NUM] [double mass] [int potential number] [double tol] [double res]


################################################################################

# Insert the desired parameters for your simulation
# If you want more than one parameter use multiple parameters by ('' '' '')
# For different variation than NUM the SIZE has to be changed
NUM_VAR=('13')
M_VAR=('2.35')
POT='0'
TOL_VAR=('1e-10')
RES_VAR=('1e-10')

# chane SIZE if necessary by commenting in/out
SIZE=${NUM_VAR[@]}
# SIZE=${M_VAR[@]}
# SIZE=${TOL_VAR[@]}
# SIZE=${RES_VAR[@]}

for (( i=0 ; i < $SIZE ; i++ )); do
  NUM=${NUM_VAR[i}
  M=${M_VAR[i}
  TOL=${TOL_VAR[i}
  RES=${RES_VAR[i}
  echo "RUNNING your simulation number: ${i}!!"
  ../main.o ${NUM} ${M} ${POT} ${TOL} ${RES}
