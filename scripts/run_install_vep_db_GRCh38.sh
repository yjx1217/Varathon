
source ./../env.sh

vep_db_dir="$VARATHON_HOME/data/vep_db/"
species_tag="homo_sapiens"
assembly_tag="GRCh38"
$vep_dir/vep_install -a cf --SPECIES $species_tag -ASSEMBLY $asembly_tag -c $vep_db_dir/$species_tag/$assembly_tag/ --CONVERT
