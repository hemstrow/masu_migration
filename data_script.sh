#!/bin/bash -l
#SBATCH -t 4-12:00:00
#SBATCH --mem 20G

printf 'Sample\t\t\t\t\t\t\t\tMapped\t\tProperly Paired\t\tSingletons\n' #print the first line: Sample  Mapped Properly Paired Singletons \n means next line
for f in *.samtools #loop for all files(f) ending with .samtools
do
    {
	printf "%s\t\t" "$f"
	let formatter=1 #acts as a counter for the first percentage value (%mapped)
	for i in `grep -E -o '\<[0-9]{1,2}\.[0-9]{2,5}\>%' $f`; do #grep looks for the decimals and i=the percentage
	    printf "%s\t\t" "$i"
            if [[ $formatter = 2 ]] #counter for the second percentage 
	        then
  		    printf "\t" #print the percentage + tab(space)
	    fi
	    ((formatter++)) #this if statement looks for the 2nd % and recognizes it because of the extra spaces after the 2nd %
 	done
	printf "\n" #print the next line(go to next line)
    }
done

