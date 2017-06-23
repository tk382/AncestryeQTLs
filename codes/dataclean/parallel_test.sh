module restore
module load gcc/6.2.0
module load parallel
cd /group/im-lab/nas40t2/tae/differentialeQTLs/codes/dataclean
filename="/group/im-lab/nas40t2/tae/differentialeQTLs/data/tissues"
echo $filename
echo -r $filename
parallel "name={}; echo $name; Rscript cleanDataForLamp.R $name" ::: line

#while read -r line
#do
#        max_bg_procs 6
#        name="$line"
#        Rscript cleanDataForLamp.R $name &
#done < "$filename"



