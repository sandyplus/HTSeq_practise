import HTSeq
fastq_file = HTSeq.FastqReader( "yeast_RNASeq_excerpt_sequence.txt", "solexa" )


fastq_file


import itertools
for read in itertools.islice( fastq_file, 10 ):
    print read


read


read.name

read.seq

read.qual



import numpy
len( read )

qualsum = numpy.zeros( len(read), numpy.int )


nreads = 0
for read in fastq_file:
    qualsum += read.qual
    nreads += 1

qualsum / float(nreads)



from matplotlib import pyplot
pyplot.plot( qualsum / nreads )    

pyplot.show()                      


alignment_file = HTSeq.SAM_Reader( "yeast_RNASeq_excerpt.sam" )


nreads = 0
for aln in alignment_file:
    qualsum += aln.read.qual
    nreads += 1


bam_reader = HTSeq.BAM_Reader( "SRR001432_head_sorted.bam" )
for a in itertools.islice( bam_reader, 5 ):  # printing first 5 reads
    print a


bam_writer = HTSeq.BAM_Writer.from_BAM_Reader( "region.bam", bam_reader ) #set-up BAM_Writer with same header as reader
for a in bam_reader.fetch( region = "1:249000000-249200000" ): #fetching reads in a region
    print "Writing Alignment", a, "to file", bam_writer.filename
    bam_writer.write( a )    


bam_writer.close()


aln

aln.read
aln.read.name
aln.read.seq
aln.read.qual


aln.iv
aln.iv.chrom
aln.iv.start
aln.iv.end
aln.iv.strand


chromlens = { 'chr1': 3000, 'chr2': 2000, 'chr1': 1000 }
ga = HTSeq.GenomicArray( chromlens, stranded=False, typecode="i" )
iv = HTSeq.GenomicInterval( "chr1", 100, 120, "." )
ga[iv] = 5
iv = HTSeq.GenomicInterval( "chr1", 110, 135, "." )
ga[iv] += 3
iv = HTSeq.GenomicInterval( "chr1", 90, 140, "." )
list( ga[iv] )  



for iv2, value in ga[iv].steps():
    print iv2, value

cvg = HTSeq.GenomicArray( "auto", stranded=True, typecode="i" )

alignment_file = HTSeq.SAM_Reader( "yeast_RNASeq_excerpt.sam" )
cvg = HTSeq.GenomicArray( "auto", stranded=True, typecode='i' )
for alngt in alignment_file:
    if alngt.aligned:
       cvg[ alngt.iv ] += 1


pyplot.plot( list( cvg[ HTSeq.GenomicInterval( "III", 200000, 500000, "+" ) ] ) )     


cvg.write_bedgraph_file( "plus.wig", "+" )
cvg.write_bedgraph_file( "minus.wig", "-" )


gas = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
gas[ HTSeq.GenomicInterval( "chr1", 100, 250 ) ] += "A"
gas[ HTSeq.GenomicInterval( "chr1", 360, 640 ) ] += "A"
gas[ HTSeq.GenomicInterval( "chr1", 510, 950 ) ] += "B"


read_iv = HTSeq.GenomicInterval( "chr1", 450, 800 )


for iv, val in gas[ read_iv ].steps():
    print iv, val



fset = set()
for iv, val in gas[ read_iv ].steps():
    fset = val

print fset



reduce( set.union, ( val for iv, val in gas[ read_iv ].steps() ) )


gtf_file = HTSeq.GFF_Reader( "Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz",
    end_included=True )


for feature in itertools.islice( gtf_file, 10 ):
    print feature



dir( feature )   




feature.iv
feature.source
feature.type
feature.score
feature.attr    


feature.name


exons = HTSeq.GenomicArrayOfSets( "auto", stranded=False )


for feature in gtf_file:
    if feature.type == "exon":
       exons[ feature.iv ] += feature.name




iv = HTSeq.GenomicInterval( "III", 23850, 23950, "." )


list( exons[iv].steps() )   






iset = None
for iv2, step_set in exons[iv].steps():
    if iset is None:
       iset = step_set.copy()
    else:
       iset.intersection_update( step_set )

print iset



list(iset)[0]



counts = {}
for feature in gtf_file:
    if feature.type == "exon":
       counts[ feature.name ] = 0


sam_file = HTSeq.SAM_Reader( "yeast_RNASeq_excerpt.sam" )
for alnmt in sam_file:
    if alnmt.aligned:
       iset = None
       for iv2, step_set in exons[ alnmt.iv ].steps():
           if iset is None:
              iset = step_set.copy()
           else:
              iset.intersection_update( step_set )
       if len( iset ) == 1:
          counts[ list(iset)[0] ] += 1


for name in sorted( counts.keys() ):
    print name, counts[name]   


















