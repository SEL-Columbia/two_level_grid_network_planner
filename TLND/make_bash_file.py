import os

ordered_files = sorted(os.listdir('../../iterative_split/'))
subfile= []
for root, dirs, files in os.walk('../../iterative_split/'):
    for name in files:
        if 'MV' not in root and 'FinalGrid' not in root:
            if 'scale' in root and '.shp' in name:
                subfile.append(os.path.join(root,name))



output =[]
for i in range(len(ordered_files)):
    item = ordered_files[i]
    cur_ward =[j[:-4] for j in subfile if item in j]
    for w in cur_ward:
        subgrid = int((w.split('/')[-1]).split('_')[1])
        further_split = (w.split('/')[-1]).split('_')[-1]
        output.append([i, subgrid,further_split])


f= open("sample_sbatch.sh","w")
f.write("#!/bin/bash\n")
f.write("\n")
for item in output:
    f.write("sbatch --mem=3500 sjob-tlnd")
    f.write(" ")
    f.write(str(item[0]))
    f.write(" ")
    f.write(str(item[1]))
    f.write(" ")
    f.write(item[2])
    f.write("\n")
f.close()

