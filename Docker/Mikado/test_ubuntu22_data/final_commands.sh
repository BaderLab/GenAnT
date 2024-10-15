#commands used that worked - with github data


echo "Mikado configure"
#configure
docker run -it -v  "$(pwd)":/global baderlab/mikado:ubuntu22_mikado2.3.2 mikado configure --list list.txt --reference chr5.fas --mode permissive --scoring plant.yaml  --copy-scoring plant.yaml --junctions junctions.bed -bt uniprot_sprot_plants.fasta -y configuration.yaml 

 echo "################  DONE configure   ##################"
  echo "Mikado Prepare"

#prepare
docker run -it -v  "$(pwd)":/global baderlab/mikado:ubuntu22_mikado2.3.2 mikado prepare --json-conf /global/configuration.yaml --start-method spawn -p 20 

echo "################  DONE Prepare   ##################"
 echo " "
 
 echo "Run Blast"

#run blast  
 docker run -it -v  "$(pwd)":/global  ncbi/blast:2.11.0 makeblastdb -in /global/uniprot_sprot_plants.fasta -dbtype prot -parse_seqids > blast_prepare.log
 
 docker run -it -v  "$(pwd)":/global  ncbi/blast:2.11.0  blastx -max_target_seqs 5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop" -num_threads 10 -query /global/mikado_prepared.fasta -db /global/uniprot_sprot_plants.fasta -out /global/mikado_prepared.blast.tsv
 
 echo "################  DONE Blast  ##################"
 echo " "
 
 echo "Run Prodigal"

 #run prodigal 
 docker run -it -v  "$(pwd)":/global metashot/prodigal:2.6.3-1 prodigal -i /global/mikado_prepared.fasta -g 1 -o /global/mikado.orfs.gff3 -f gff
 
 echo "################  DONE Prodigal   ##################"
 echo " "
  
 echo "Run Serialize"

 
 #serialize
 #not working with - --tsv /global/mikado_prepared.blast.tsv --orfs mikado.orfs.gff3  -bt /global/uniprot_sprot_plants.fasta
docker run -it -v  "$(pwd)":/global baderlab/mikado:ubuntu22_mikado2.3.2 mikado serialise --json-conf /global/configuration.yaml  --junctions junctions.bed 

 echo "################  DONE Serialise   ##################"
 echo " "
  
 echo "Run Pick"


#pick
docker run -it -v  "$(pwd)":/global baderlab/mikado:ubuntu22_mikado2.3.2 mikado pick --json-conf /global/configuration.yaml -db /global/mikado.db  --mode permissive /global/mikado_prepared.gtf --loci-out /global/pick_annotation.out --log /global/mikado_pick.log

 echo "################  DONE Pick  ##################"
 echo " "
  

