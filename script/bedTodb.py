#! /bin/env python

import sqlite3
import sys


db = sys.argv[2]
conn = sqlite3.connect(db)
c = conn.cursor()

#open bedfile
bedfile = open(sys.argv[1],'r')

c.execute('CREATE TABLE chromosomes(id text, length int)')

# add each bedfile line to db
for line in bedfile:
    [chrom, start, end, info] = line.rstrip().split('\t')
    [ID,length,span] = info.split(',')
    try: ## add line
        c.execute("INSERT INTO "+chrom+" VALUES ("+start+","+end+","+length[:-2]+",'"+ID+"')")
    except sqlite3.OperationalError: ## Create table
        c.execute("INSERT INTO chromosomes VALUES ('"+chrom+"',0)")
        c.execute('CREATE TABLE '+chrom+'(start int, end int, length int, id text)')
        c.execute("INSERT INTO "+chrom+" VALUES ("+start+","+end+","+length[:-2]+",'"+ID+"')")
bedfile.close

## Save (commit) the changes
conn.commit()

# test to fetch
for chrom in ['chrX','chrY']:
    for row in c.execute('SELECT * FROM '+chrom+' WHERE '+chrom+'.length >= 10 ORDER BY length, start'):
        print chrom, row
    
#close the connection
conn.close()
