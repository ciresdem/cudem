#!/bin/bash

function help () {
echo "create_datalist - Creates a datalist for all xyz files in a directory"
	echo "Usage: $0 datalist_name"
	echo "* datalist_name: <name of output datalist>"
}

#see if 1 parameters was provided
#show help if not
if [ ${#@} == 1 ]; 
then

	#User inputs    	
	datalist_name=$1
	ls *.xyz | awk '{print $1, 168}' > $datalist_name".datalist"
	mbdatalist -F-1 -I$datalist_name".datalist" -O -V
	#
	echo
	echo "All done"

else
	help
fi

### End




