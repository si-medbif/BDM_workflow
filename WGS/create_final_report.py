#! /usr/bin/env python

import sys
import datetime

reportfile = 'WGS_report.md'

header = '# WGS Analysis Report\n{date:%Y-%m-%d %H:%M:%S}\n'.format(date=datetime.datetime.now())

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
Figure 3: The statistics of WGS mapping results'

with open(reportfile, 'w') as fout:
    fout.write(header)
    fout.write('\n\n')
    fout.write(section1)
    fout.write('\n\n')
    fout.write(section2)
    fout.write('\n\n')
    fout.write(section3)
