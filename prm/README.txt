PRM solver:		src/lib/prm
Geo lib:		src/lib/collision4
Robot application:	src/lib/robotspace
Array configspace:	src/lib/arrayspace
application with array:	src/array.cpp
application with robot:	src/main.cpp


Using the robot application:

- generate robot configuration with matlab/saveconfig.m
- run ./main or ./maincuda (with or without GPU)
- plot the result with matlab/robotplot.m

Using the array application:

- store int array in binary file: array.bin (already done)
- run ./array
- plot with matlab/graphplot.m

