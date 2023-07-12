# Some COVID Analysis
Calculates the mean and standard deviation of the lengths and gc content of coronavirus genomes
Coronavirus geneomes downloaded from NCBI (see `data/data.md` for more details)

**1. Load in Data**
```julia
using BISC195
genomes = parse_fasta("../data/cov-sequences.fasta")
``` 

**2. Calculate the mean length of the genomes**
```julia
len = 0
count = length(genomes[2])
for sequence in genomes[2]
    len += length(sequence)
end
mean_len = len/count
```

**3. Calculate the standard deviation of the lengths of the genomes**
```julia
distance = 0
for sequence in genomes[2]
    distance += (length(sequence)-mean_len)^2
end
std_len = sqrt(distance/(count-1))
```

**4. Calculate the gc content of the genomes**
```julia
gcs = []
for sequence in genomes[2]
   push!(gcs, gc_content(sequence))
end
```

**5. Calculate the mean gc content of the genomes**
```julia
total = 0
for gc in gcs
    total += gc
end
mean_gc = total/count
```

**6. Calculate the the standard deviation of the gc content of the genomes**
```julia
distance = 0
for gc in gcs
    distance += (gc-mean_gc)^2
end
std_gc = sqrt(distance/(count-1))
```