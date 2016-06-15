#!/bin/bash -e

# curl -O ftp://hgdownload.soe.ucsc.edu//apache/htdocs/goldenPath/hg19/database/snp137Common.txt.gz
# curl -O ftp://hgdownload.soe.ucsc.edu//apache/htdocs/goldenPath/hg19/database/snp135Common.txt.gz
# curl -O ftp://hgdownload.soe.ucsc.edu//apache/htdocs/goldenPath/hg19/database/snp132Common.txt.gz

gunzip -c snp137Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp137Common_small.txt.gz
gunzip -c snp135Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp135Common_small.txt.gz
gunzip -c snp132Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp132Common_small.txt.gz

