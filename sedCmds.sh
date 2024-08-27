#! /usr/bin/bash
time sed -n '1~4s/^@/>/p;2~4p' input10MB.fastq > sedFASTA10MB.fasta
time sed -n '1~4s/^@/>/p;2~4p' input100MB.fastq > sedFASTA100MB.fasta
time sed -n '1~4s/^@/>/p;2~4p' input1GB.fastq > sedFASTA1GB.fasta
time sed -n '1~4s/^@/>/p;2~4p' input8GB.fastq > sedFASTA8GB.fasta