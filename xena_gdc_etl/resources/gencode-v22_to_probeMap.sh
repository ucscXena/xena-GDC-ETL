#!/bin/bash

echo -e "id\tgene\tchrom\tchromStart\tchromEnd\tstrand" > gencode.v22.annotation.gene.probeMap
awk 'BEGIN { FS="\t"; OFS="\t" } $3 == "gene" { match($9, /gene_id\s\"([^\"]+)\";/, gid); match($9, /gene_name\s\"([^\"]+)\";/, gname); print gid[1], gname[1], $1, $4, $5, $7 }' gencode.v22.annotation.gtf >> gencode.v22.annotation.gene.probeMap
