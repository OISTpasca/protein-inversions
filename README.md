# protein-inversions
Pipeline to detect protein residue inversions in a duplication phylogeny

## scoring function
We start from four aa-frequency arrays [a b c d] of length 20. The four arrays represent a partition of the
 multiple sequence alignment in two groupings of two. One grouping represent the gene duplication event. 
 [a b] group sequences of protein duplicate 1, while [c d] group sequences of protein duplicate 2.
The second grouping represent a subset of species. [a d] group different species than [b c].  
The probability of two similar events is being calculated. **META** is the event of an inversion between the duplicates 
in a subgroup of species. So to say, that [b c] are inverted and the highest similarity is across duplicate proteins.
**SS**, for species specific (adaptation) is the event of highest similarity among the same species, or a virtual 
inversion between [a b].  
To calculate the probability of these events we proceed as follows:  
We define the probability of an amino acid match between two groupings as
```math
P_{(a=b)} = \sum_{i}^{20} a_{i} b_{i}
```
We assume that the probability of matching of two arrays is independent:
```math
P_{(a=b,c=d)} = P_{(a=b)} * P_{(c=d)}
```
Then, we calculate the conditional probability of a match on a different amino acid given the double match:
```math
P_{(b!=c|a=c,b=d)} = \sum_{i}^{20} \frac{a_i * c_i}{P_{(a=c)}}* (1-\frac{b_i * d_i}{P_{(b=d)}})
```
for the **META** event and:
```math
P_{(b!=c|a=d,b=c)} = \sum_{i}^{20} \frac{a_i * d_i}{P_{(a=d)}}* (1-\frac{b_i * c_i}{P_{(b=c)}})
```
for the **SS** event.  
To obtain the joined probability, we just multiply by the probability of the condition. For example, for **META**:
```math
P_{(b!=c,a=c,b=d)} = P_{(b!=c|a=c,b=d)} * P_{(a=c,b=d)}
```
