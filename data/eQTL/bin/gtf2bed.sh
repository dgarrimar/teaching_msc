cat $1 | perl -ne 'my @l=split(/\t/,$_);if($l[2] eq "gene"){$_ =~/gene_id "([^"]+)"/;print $l[0]."\t".($l[3]-1)."\t".$l[4]."\t".$1."\t".$l[5]."\t".$l[6]."\n";}'
