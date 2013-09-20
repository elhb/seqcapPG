#! /bin/env python

import sqlite3
import sys
import pysam

color = True
if color:bs = "\033[1m";be = "\033[0;0m";PURPLE = '\033[95m';BLUE = '\033[94m';GREEN = '\033[92m';YELLOW = '\033[93m';RED = '\033[91m';ENDC = '\033[0m';BLACK = '\33['+'0;30'+'m';CYAN = '\33['+'0;36'+'m';GRAY = '\33['+'0;37'+'m';BBLUE = '\33['+'1;34'+'m';BRED = '\33['+'1;31'+'m';BGREEN = '\33['+'1;32'+'m';BCYAN = '\33['+'1;36'+'m';BPURPLE = '\33['+'1;35'+'m';BYELLOW = '\33['+'1;33'+'m';DGRAY = '\33['+'1;30'+'m'
else:bs = '';be = '';HEADER = '';BLUE = '';GREEN = '';YELLOW = '';RED = '';ENDC = '';PURPLE='';BLACK='';CYAN='';GRAY='';BBLUE='';BRED='';BGREEN='';BCYAN='';BPURPLE='';BYELLOW='';DGRAY='';

plusminus = 100#bp
endge_size = 5
bamfiles = [
    'TEMPO/marta/marta.i1.c15.recal.final.bam',
    'TEMPO/marta/marta.i2.c21.recal.final.bam',
    'TEMPO/marta/marta.i3.c22.recal.final.bam',
    'TEMPO/marta/marta.i4.c24.recal.final.bam',
    'TEMPO/marta/marta.i5.c29.recal.final.bam',
    'TEMPO/marta/marta.i6.c30.recal.final.bam',
    'TEMPO/marta/marta.i7.c38.recal.final.bam',
    'TEMPO/marta/marta.i8.c40.recal.final.bam',
    'TEMPO/marta/marta.i9.genomic.recal.final.bam',
]

def main():
    referencefa=pysam.Fastafile(sys.argv[3])
    db = sys.argv[2]
    conn = sqlite3.connect(db)
    c = conn.cursor()
    total_overlaps = 0
    for chrom in c.execute('SELECT id FROM chromosomes').fetchall(): #[['chrY','']]:#
        chrom = str(chrom[0])
        print chrom
        for row in c.execute('SELECT * FROM '+chrom+' WHERE '+chrom+'.length >= 5 ORDER BY length, start'):
            [start, end,length,ID] = row
#            if start < 43421992-100 or end > 43421992+100: continue
            print '\n###############################\n'+ID+', '+chrom+':'+str(start)+'-'+str(end)+', '+str(length)+'bp'
            for bamname in bamfiles:
                overlaps = processbam(row,bamname,chrom,referencefa)#sys.argv[1]
                total_overlaps += overlaps
    #bam.close()
    conn.close()
    referencefa.close()
    print 'we found', total_overlaps, 'overlapping reads'

def processbam(row,bamname,chrom,referencefa):
        [start, end,length,ID] = row
        start=start-1
        ID = str(ID)

        #get count of reads overlapping region
        bam = pysam.Samfile(bamname,'rb')
        header = '' 
        header += bamname +' '+str(bam.count(chrom, start-plusminus-1, end+plusminus))+' reads within +-'+str(plusminus)+'bp,'

        overlaps = 0
        print header

        # initiate inforamtion holders
        positionsArray = []
        by_relpos = {}

        #Fill in reference data and relative position names
        ref=referencefa.fetch(chrom, start-endge_size, end+endge_size+1)
        for i in range(len(ref)):
            current_pos = i+start-endge_size
            by_relpos[current_pos-start] = {'Reference_Sequence':ref[i]}
            if current_pos-start <= length-1:   relpos_name = str(current_pos-start)
            else:                               relpos_name = '+'+str(current_pos-start-length+1)
            while len(relpos_name) < 2: relpos_name+=' '
            by_relpos[current_pos-start]['Relative_position']= relpos_name

        # Fill in per base information
        indel = False
        for col in bam.pileup(chrom, start-endge_size-1, end+endge_size-1):
            if col.pos >= start-endge_size and col.pos <= end+endge_size:
                current_column = {'position':col.pos,'readDepth':col.n}
                for pileupread in col.pileups:
                    
                    relpos=col.pos-start
                    current_column['relpos'] = relpos
                    
                    if   pileupread.alignment.is_read1:read_number = '_1'
                    elif pileupread.alignment.is_read2:read_number = '_2'
                    else: print'ERROR-SCHMERROR!!!!',sys.exit
                    
                    by_relpos[relpos][pileupread.alignment.qname+read_number] = pileupread.alignment.seq[pileupread.qpos]
                    current_column[pileupread.alignment.qname+read_number] = pileupread.alignment.seq[pileupread.qpos]
                    
                    if pileupread.is_del:
                        current_column[pileupread.alignment.qname+read_number] = '-'
                        by_relpos[relpos][pileupread.alignment.qname+read_number] = '-'
                    #if pileupread.indel and relpos>= 0 and relpos <= length+1: indel = True #print pileupread.indel, 'BASE INDEL!!! ---- at ',col.pos-start,'in', pileupread.alignment.qname, pileupread.alignment.seq
                    if pileupread.indel > 0:
                        #current_column[pileupread.alignment.qname] = str(pileupread.indel)
                        #by_relpos[relpos][pileupread.alignment.qname] = str(pileupread.indel)
                        by_relpos[relpos][pileupread.alignment.qname+read_number] = RED+''.join([pileupread.alignment.seq[i] for i in range(pileupread.qpos,pileupread.qpos+pileupread.indel+1,1)])+ENDC
                        current_column[pileupread.alignment.qname+read_number] = ''.join([pileupread.alignment.seq[i] for i in range(pileupread.qpos,pileupread.qpos+pileupread.indel+1,1)])
                positionsArray.append(current_column)

        #get all read names
        reads = {}
        for position in positionsArray:
            for read_name in position.keys():
                if read_name != 'readDepth' and read_name != 'position' and read_name != 'relpos': reads[read_name] = True

        positions = by_relpos.keys()
        positions.sort()

        # delete reads with lack of coverage and/or missmatches in the "outer edge regions"
        for read_name in reads.keys():
                for position in range(-1*endge_size,-1)+range(length+1,length+endge_size):
#                    print position
                    try: out = by_relpos[position][read_name]
                    except KeyError:
                        del reads[read_name];
                        break
                    if out != by_relpos[position]['Reference_Sequence']:
                        del reads[read_name];
                        break

        #get consensus
        for position in positionsArray:
            by_relpos[position['relpos']]['ConsensusSequence'] = ''
            base_count = {'total':0}

            for read_name in reads.keys():
                if read_name != 'readDepth' and read_name != 'position' and read_name != 'relpos':
                    base = position[read_name]
                    base_count['total']+=1
                    try: base_count[base] +=1
                    except KeyError: base_count[base] =1


            for base,count in base_count.iteritems():
                if base == 'total': continue
                by_relpos[position['relpos']]['ConsensusSequence'] += base+'='+str(round(100*float(count)/base_count['total'],0))+','#+'%,'
            #by_relpos[position['relpos']]['ConsensusSequence'] += ' '


        # Print graphical output
        for read_name in ['Relative_position','Reference_Sequence']+reads.keys()+['ConsensusSequence']:
        #for read_name in ['Relative_position','Reference_Sequence']+['ConsensusSequence']:
                print read_name+'\t',
                if read_name =='Relative_position': print '\t\t\t',
                if read_name =='Reference_Sequence':print '\t\t\t',
                if read_name =='ConsensusSequence':
                    print '\t\t\t',
                    for position in positions:
                        try :bases = by_relpos[position][read_name]
                        except KeyError: bases = None
                        if bases:
                            bases = bases.split(',')
                            tmp = {}
                            out = ''
                            #print 'bases->',bases
                            for base in bases[:-1]:
                                base = base.split('=')
                                perc = float(base[1])
                                base = base[0]
                                tmp[base] = perc
                            for base,perc in tmp.iteritems():
                                if perc > 65.0: out = base
                                if perc > 35.0 and perc < 65.0:
                                    if out: out += '/'+base
                                    else:out += base
                        else:out = 'NA'

                        if out != by_relpos[position]['Reference_Sequence'] and read_name != 'Relative_position' and out != '.' and out != 'NA': color = RED
                        else: color =''
                        if len(out) < 2: out +=' '
                        print color+out+ENDC+'\t',
                    print ''                    
                    continue

                for position in positions:

                    try: out = by_relpos[position][read_name]
                    except KeyError: out = '.'

                    if out != by_relpos[position]['Reference_Sequence'] and read_name != 'Relative_position' and out != '.': color = RED
                    else: color =''

                    if len(out) < 2: out +=' '
                    
                    print color+out+ENDC+'\t',

                print ''

        return overlaps
    
#####
#check if run or imported // call main() or not
#####
if __name__ == "__main__":
    main()
#END of script