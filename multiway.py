import sys
import pandas as pd

infile=sys.argv[1]
out=sys.argv[2]
outname=sys.argv[3]

dat = pd.read_csv(infile, sep="\t", header=None)

dat.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'strand1', 'strand2', 'placeholder', 'readid','chrom_probe','start_probe', 'end_probe', 'probeID', 'probeName']

##remove duplicated rows with the same mapping location, readid and probe regions
newdf = dat.drop_duplicates(
  subset = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'readid','chrom_probe', 'start_probe', 'end_probe'],
  keep = 'first').reset_index(drop = True)

##keep the rows with duplicated readid but with different mapping regions and probeID. 
duplicate = newdf[newdf.duplicated(["readid"], keep=False)]

duplicate.to_csv(f'{outname}_multiway.txt', sep='\t', index=False)

### remove duplicated rows with the same readid and probeName.
newdf2 = dat.drop_duplicates(subset=['readid', 'probeName'], keep='first').reset_index(drop=True)
duplicate = newdf2[newdf2.duplicated(["readid"], keep=False)]
### output probe combinations for the same readid.
probeComb = duplicate.groupby(['readid']).agg({'probeName':','.join}).reset_index()
probeComb.to_csv(f'{outname}_probeCombination.txt', sep='\t', index=False)

#count = {}
#for index, row in probeComb.iterrows():
#    probeSet = tuple(row['probeName'].split(","))
#    
#    sortedSet = tuple(sorted(probeSet))
#    #print(sortedSet)
#    
#    if sortedSet not in count.keys():
#        count[sortedSet] = 1
#    else:
#        count[sortedSet] += 1

## filter for multiway count 
#for i in count.copy():
#    if len(i)<=2:
#        del count[i]
        

count={}
twoway_read=[]
threeway_read=[]
fourway_read=[]
multiway_read=[]

for index, row in probeComb.iterrows():
    probeSet = set(tuple(row['probeName'].split(",")))
    
    sortedSet = tuple(sorted(probeSet))
    #print(sortedSet)
    
    if sortedSet not in count.keys() :
        count[sortedSet] = 1
        if len(sortedSet) == 2:
            twoway_read.append(row['readid'])
        elif len(sortedSet) == 3:
            threeway_read.append(row['readid'])
            multiway_read.append(row['readid'])
        elif len(sortedSet) > 3:
            fourway_read.append(row['readid'])
            multiway_read.append(row['readid'])
        else:
            pass
    else:
        count[sortedSet] += 1
        if len(sortedSet) == 2:
            twoway_read.append(row['readid'])
        elif len(sortedSet) == 3:
            threeway_read.append(row['readid'])
            multiway_read.append(row['readid'])
        elif len(sortedSet) > 3:
            fourway_read.append(row['readid'])
            multiway_read.append(row['readid'])
        else:
            pass
    #break

multiway_interaction = dat[dat['readid'].isin(multiway_read)]
multiway_interaction.to_csv(f'{outname}_multiread.info.bed', sep="\t", index=False)

twoway_interaction = dat[dat['readid'].isin(twoway_read)]
twoway_interaction.to_csv(f'{outname}_twoway_read.info.bed', sep="\t", index=False)

#threeway_interaction = dat[dat['readid'].isin(twoway_read)]
#threeway_interaction.to_csv(f'{outname}_threeway_read.info.bed', sep="\t", index=False)

fourway_interaction = dat[dat['readid'].isin(fourway_read)]
fourway_interaction.to_csv(f'{outname}_fourway_read.info.bed', sep="\t", index=False)


multiComb = multiway_interaction.groupby(['readid']).agg({'probeName':','.join}).reset_index()
multiComb.to_csv(f'{outname}_multiCombination.txt', sep='\t', index=False)

twowayComb = twoway_interaction.groupby(['readid']).agg({'probeName':','.join}).reset_index()
twowayComb.to_csv(f'{outname}_twowayCombination.txt', sep='\t', index=False)

#threewayComb = twoway_interaction.groupby(['readid']).agg({'probeName':','.join}).reset_index()
#threewayComb.to_csv(f'{outname}_threewayCombination.txt', sep='\t', index=False)

fourwayComb = fourway_interaction.groupby(['readid']).agg({'probeName':','.join}).reset_index()
fourwayComb.to_csv(f'{outname}_fourwayCombination.txt', sep='\t', index=False)


multi_count = {}
for index, row in multiComb.iterrows():
    probeSet = set(row['probeName'].split(","))
    
    sortedSet = tuple(sorted(probeSet))
    #print(sortedSet)
    
    if sortedSet not in multi_count.keys():
        multi_count[sortedSet] = 1
    else:
        multi_count[sortedSet] += 1
        
stat = pd.DataFrame(multi_count.items())
stat.columns=['combination', 'number_of_reads']
stat.to_csv(f'{outname}_multi_stats.txt', sep='\t', index=False)

multi =multiway_interaction.drop_duplicates(subset = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'readid'],keep = 'first').reset_index(drop = True)
multiDup = multi[multi.duplicated(["readid"], keep=False)]
#multiDup.to_csv(f'{outname}_multiread.info.bed', sep="\t", index=False) #no need to output this. 


multiBed = []
n=0
readdict = {}

for index, row in multiDup.iterrows():
    
    if row['readid'] not in readdict.keys():
        readdict[row['readid']] = n
        multiBed.append([row['chrom1'], row['start1'], row['end1'], row['readid'], n, "."])
        multiBed.append([row['chrom2'], row['start2'], row['end2'], row['readid'], n, "."])
        n +=1
    else:
     #   readdict[row['readid']] += 1
        multiBed.append([row['chrom1'], row['start1'], row['end1'], row['readid'], readdict[row['readid']], "."])
        multiBed.append([row['chrom2'], row['start2'], row['end2'], row['readid'], readdict[row['readid']], "."])

multiBed_DF = pd.DataFrame(multiBed)
multiBed_DF.columns=['chrom','start','end','readid','group', 'strand']
multiBed_DF = multiBed_DF.drop_duplicates(
  subset = ['chrom', 'start', 'end', 'readid', 'group', 'strand'],
  keep = 'first').reset_index(drop = True)

multiBed_DF.to_csv(f'{outname}_multiread.forviz.bed', sep="\t", index=False, header=False)

#stat = pd.DataFrame(count.items())
#stat.columns=['combination', 'number_of_reads']
#stat.to_csv(f'{outname}_stats.txt', sep='\t', index=False)


#### output four-way interactions for visulization
multi =fourway_interaction.drop_duplicates(subset = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'readid'],keep = 'first').reset_index(drop = True)
multiDup = multi[multi.duplicated(["readid"], keep=False)]
#multiDup.to_csv(f'{outname}_multiread.info.bed', sep="\t", index=False) #no need to output this. 


multiBed = []
n=0
readdict = {}

for index, row in multiDup.iterrows():
    
    if row['readid'] not in readdict.keys():
        readdict[row['readid']] = n
        multiBed.append([row['chrom1'], row['start1'], row['end1'], row['readid'], n, "."])
        multiBed.append([row['chrom2'], row['start2'], row['end2'], row['readid'], n, "."])
        n +=1
    else:
     #   readdict[row['readid']] += 1
        multiBed.append([row['chrom1'], row['start1'], row['end1'], row['readid'], readdict[row['readid']], "."])
        multiBed.append([row['chrom2'], row['start2'], row['end2'], row['readid'], readdict[row['readid']], "."])

multiBed_DF = pd.DataFrame(multiBed)
multiBed_DF.columns=['chrom','start','end','readid','group', 'strand']
multiBed_DF = multiBed_DF.drop_duplicates(
  subset = ['chrom', 'start', 'end', 'readid', 'group', 'strand'],
  keep = 'first').reset_index(drop = True)

multiBed_DF.to_csv(f'{outname}_fourwayread.forviz.bed', sep="\t", index=False, header=False)