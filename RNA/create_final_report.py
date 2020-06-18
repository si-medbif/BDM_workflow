#! /usr/bin/env python

import sys
import datetime

reportfile = 'RNA_sequencing_report.md'

header = '# RNA-Seq Analysis Report\n{date:%Y-%m-%d %H:%M:%S}\n'.format(date=datetime.datetime.now())

section1 = '## 1. Summary\n \
This report includes quality check for raw sequencing data (section 2), \
reads mapping and assignment (section 3).'

section2 = '## 2. Sequencing Quality\n \
### 2.1 Summary of sequencing quality\n \
![Base quality](figure1.png) \n \
Figure 1: The mean quality scores of sequencing reads in each position \n \
![GC scores](figure2.png) \n \
Figure 2: The average GC content of sequencing reads'

section3 = '## 3. Mapping quality\n \
### 3.1 Summary of mapping quality\n \
![Reads mapping](figure3.png) \n \
Figure 3: The statistics of RNAseq mapping results'

df = pd.read_table('table1.tsv', header=0)
df['Ratio'] = df['Ratio'].round(decimals=3)
table1 = 'Table1: The summary of RNAseq mapping results\n\n|'
section = '|'
for c in df.columns:
    table1 += c+'|'
    section += '---'+'|'
table1 += '\n' + section + '\n'
for ind, row in df.iterrows():
    line = '|'
    for c in row:
        line += str(c)+'|'
    table1 += line + '\n'

with open(reportfile, 'w') as fout:
    fout.write(header)
    fout.write('\n\n')
    fout.write(section1)
    fout.write('\n\n')
    fout.write(section2)
    fout.write('\n\n')
    fout.write(section3)
    fout.write('\n\n')                                                                                                                                                       fout.write(table1)   
