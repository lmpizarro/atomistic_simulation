#!/bin/bash

array=(masa_1.0 masa_2.0 masa_3.0 masa_4.0 masa_5.0 masa_6.0 masa_7.0 masa_8.0 masa_0.6 masa_1.6 masa_2.6 masa_3.6 masa_4.6 masa_5.6 masa_6.6 masa_7.6)


echo "Array size: ${#array[*]}"

echo "Array items:"
for item in ${array[*]}
do
    printf "   %s\n" $item
    cd $item
    ls
    qsub run.qsub
    cd ..
done
