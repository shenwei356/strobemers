# Evaluation

## Number of matched strobemers

query: 150bp, snp: 1 (0.006)

    $ go run test1_matches.go  q0-snp1.fasta r0.fasta  | csvtk pretty -t
    query     ref   method                           nQuery   nRef      nCommon   qCov
    -------   ---   ------------------------------   ------   -------   -------   -----
    q0-snp1   r0    Kmer(20)                         131      1546586   111       84.73 *
    q0-snp1   r0    MinStrobes(2,10,12,12,shrink)    131      1548767   109       83.21
    q0-snp1   r0    MinStrobes(2,10,12,12)           129      1548765   109       84.50
    q0-snp1   r0    RankStrobes(2,10,12,12,shrink)   131      1548767   109       83.21
    q0-snp1   r0    RankStrobes(2,10,12,12)          129      1548765   109       84.50
                                                                                
    q0-snp1   r0    Kmer(21)                         130      1547218   109       83.85
    q0-snp1   r0    MinStrobes(3,7,9,9,shrink)       126      1549315   108       85.71 *
    q0-snp1   r0    MinStrobes(3,7,9,9)              126      1549315   108       85.71 *
    q0-snp1   r0    RankStrobes(3,7,9,9,shrink)      126      1549315   108       85.71 *
    q0-snp1   r0    RankStrobes(3,7,9,9)             126      1549315   108       85.71 *
                                                                                
    q0-snp1   r0    Kmer(20)                         131      1546586   111       84.73
    q0-snp1   r0    MinStrobes(2,10,12,16,shrink)    131      1548376   107       81.68
    q0-snp1   r0    MinStrobes(2,10,12,16)           125      1548370   107       85.60  *
    q0-snp1   r0    RankStrobes(2,10,12,16,shrink)   131      1548438   108       82.44
    q0-snp1   r0    RankStrobes(2,10,12,16)          125      1548432   108       86.40  *
                                                                                
    q0-snp1   r0    Kmer(21)                         130      1547218   109       83.85
    q0-snp1   r0    MinStrobes(3,7,9,13,shrink)      122      1545403   102       83.61
    q0-snp1   r0    MinStrobes(3,7,9,13)             118      1545399   102       86.44 **
    q0-snp1   r0    RankStrobes(3,7,9,13,shrink)     122      1545522   107       87.70 **
    q0-snp1   r0    RankStrobes(3,7,9,13)            118      1545518   105       88.98 **
    
query: 150bp, snp: 3 (0.02)

    $ go run test1_matches.go  q2-snp3.fasta r2.fasta  | csvtk pretty -t
    query     ref   method                           nQuery   nRef      nCommon   qCov
    -------   ---   ------------------------------   ------   -------   -------   -----
    q2-snp3   r2    Kmer(20)                         131      1687558   84        64.12 *
    q2-snp3   r2    MinStrobes(2,10,12,12,shrink)    131      1687785   82        62.60
    q2-snp3   r2    MinStrobes(2,10,12,12)           129      1687783   82        63.57
    q2-snp3   r2    RankStrobes(2,10,12,12,shrink)   131      1687785   82        62.60
    q2-snp3   r2    RankStrobes(2,10,12,12)          129      1687783   82        63.57
                                                                                
    q2-snp3   r2    Kmer(21)                         130      1687656   82        63.08
    q2-snp3   r2    MinStrobes(3,7,9,9,shrink)       126      1687865   84        66.67 *
    q2-snp3   r2    MinStrobes(3,7,9,9)              126      1687865   84        66.67 *
    q2-snp3   r2    RankStrobes(3,7,9,9,shrink)      126      1687865   84        66.67 *
    q2-snp3   r2    RankStrobes(3,7,9,9)             126      1687865   84        66.67 *
                                                                                
    q2-snp3   r2    Kmer(20)                         131      1687558   84        64.12 *
    q2-snp3   r2    MinStrobes(2,10,12,16,shrink)    131      1687487   76        58.02
    q2-snp3   r2    MinStrobes(2,10,12,16)           125      1687481   76        60.80
    q2-snp3   r2    RankStrobes(2,10,12,16,shrink)   131      1687529   77        58.78
    q2-snp3   r2    RankStrobes(2,10,12,16)          125      1687523   76        60.80
                                                                                
    q2-snp3   r2    Kmer(21)                         130      1687656   82        63.08 *
    q2-snp3   r2    MinStrobes(3,7,9,13,shrink)      122      1684611   74        60.66
    q2-snp3   r2    MinStrobes(3,7,9,13)             118      1684607   72        61.02
    q2-snp3   r2    RankStrobes(3,7,9,13,shrink)     122      1684661   71        58.20
    q2-snp3   r2    RankStrobes(3,7,9,13)            118      1684657   68        57.63
    
query: 150bp, snp: 7 (0.47)

    $ go run test1_matches.go  q1-snp7.rc.fasta r1.fasta  | csvtk pretty -t
    query        ref   method                           nQuery   nRef      nCommon   qCov
    ----------   ---   ------------------------------   ------   -------   -------   -----
    q1-snp7.rc   r1    Kmer(20)                         131      2802879   54        41.22 *
    q1-snp7.rc   r1    MinStrobes(2,10,12,12,shrink)    131      2804781   52        39.69
    q1-snp7.rc   r1    MinStrobes(2,10,12,12)           129      2804779   52        40.31
    q1-snp7.rc   r1    RankStrobes(2,10,12,12,shrink)   131      2804781   52        39.69
    q1-snp7.rc   r1    RankStrobes(2,10,12,12)          129      2804779   52        40.31
                                                                                    
    q1-snp7.rc   r1    Kmer(21)                         130      2804365   51        39.23 *
    q1-snp7.rc   r1    MinStrobes(3,7,9,9,shrink)       126      2806161   48        38.10
    q1-snp7.rc   r1    MinStrobes(3,7,9,9)              126      2806161   48        38.10
    q1-snp7.rc   r1    RankStrobes(3,7,9,9,shrink)      126      2806161   48        38.10
    q1-snp7.rc   r1    RankStrobes(3,7,9,9)             126      2806161   48        38.10
                                                                                    
    q1-snp7.rc   r1    Kmer(20)                         131      2802879   54        41.22 *
    q1-snp7.rc   r1    MinStrobes(2,10,12,16,shrink)    131      2803507   51        38.93
    q1-snp7.rc   r1    MinStrobes(2,10,12,16)           125      2803501   51        40.80
    q1-snp7.rc   r1    RankStrobes(2,10,12,16,shrink)   131      2803659   47        35.88
    q1-snp7.rc   r1    RankStrobes(2,10,12,16)          125      2803653   44        35.20
                                                                                    
    q1-snp7.rc   r1    Kmer(21)                         130      2804365   51        39.23 *
    q1-snp7.rc   r1    MinStrobes(3,7,9,13,shrink)      122      2797218   36        29.51
    q1-snp7.rc   r1    MinStrobes(3,7,9,13)             118      2797214   36        30.51
    q1-snp7.rc   r1    RankStrobes(3,7,9,13,shrink)     122      2797918   42        34.43
    q1-snp7.rc   r1    RankStrobes(3,7,9,13)            118      2797914   41        34.75
