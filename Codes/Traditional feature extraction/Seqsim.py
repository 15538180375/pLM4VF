def blast_qseqid(blastfile):
    qseqid = []
    seqid = []
    with open(blastfile) as f:
        for line in f.readlines():
            seqid.append(line.split('\t')[0])
    [qseqid.append(i) for i in seqid if not i in qseqid]
    return qseqid
