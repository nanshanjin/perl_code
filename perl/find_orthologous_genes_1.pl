#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use File::Copy;
use FindBin qw($RealBin);
my $bin = $FindBin::RealBin;
use lib ("$FindBin::Bin/PerlLib");
#use PBS::Queue;
#use Parallel::ForkManager;
use DBI;
use Cwd;

my $help;
my $pjname;
my $seqlist;
my $retainspeclist;
my $orthomcl;
my $workdir;
my $inflation_value_string;
my $percents_string;
GetOptions(
    "help!"      => \$help,
    "pjname=s"   => \$pjname,
    "seqlist=s"  => \$seqlist,
	"orthomcl!"  => \$orthomcl,
	"retained=s" => \$retainspeclist,
	"inflation=s"=> \$inflation_value_string,
	"percents=s" => \$percents_string,
    "workdir=s"  => \$workdir,
);

my $INFO=<<LINES;
Usage:
	perl $0 -orthomcl -seqlist list.txt -workdir ./ -pjname HH

Arguments:
	help        display this message.
	pjname      project name, default is "Test".
	seqlist     list file contains 3 cols: species name, typeand sequence file path.
                	The type should be orthomcl or search,
                	orthomcl: means its corresponding species should be used to find core-orthology cluster
                	search  : means its corresponding species should be find its orthologous genes
	retained    list file contains species names must be contained in each orthologous group
	inflation   inflation values used in mcl, multiple values could be speciefied and separated be ","
	            default is "1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0"
	percents    percents cut-offs used for determining orthologous gene cluster, default is "0.5,0.75,1"
	orthomcl    whether run orthomcl to find core-orthology clusters, it should be run at first time
	workdir     dir contains all tmp or res data, default is "./".

example:
	perl find_orthologous_genes.pl -pjname Test1 -seqlist seq.list -orthomcl -workdir outdir
	seq.list of example:
	Bd  orthomcl    /mnt/ilustre/users/xiaoyue.wang/software/orthologous_genes/test/Bd.faa
	Bj  orthomcl    /mnt/ilustre/users/xiaoyue.wang/software/orthologous_genes/test/Bj.faa
	C_31    orthomcl    /mnt/ilustre/users/xiaoyue.wang/software/orthologous_genes/test/C_31.faa
	F_31    orthomcl    ./F_31.faa
	Bd	search	/mnt/ilustre/users/xiaoyue.wang/software/orthologous_genes/test/Bd.faa

LINES

if($help){
	die $INFO;
}

$workdir ||= "./";
$workdir = File::Spec->rel2abs($workdir);
$pjname  ||= "Test";
$inflation_value_string ||= "1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0";
$percents_string ||= "0.5,0.75,1";

my $DIR_OrthoMCL_Fasta    = "$workdir/orthomcl_fasta";
my $DIR_OrthoMCL_DBmblast = "$workdir/orthomcl_mblast";
my $DIR_OrthoMCL_DBblast  = "$workdir/orthomcl_blast";
my $DIR_mclOutput         = "$workdir/mclOutput";
my $DIR_group             = "$workdir/group";
my $DIR_group_seq         = "$workdir/group_seq";
my $DIR_group_aln         = "$workdir/group_aln";
my $DIR_group_msf         = "$workdir/group_msf";
my $DIR_hmmer             = "$workdir/hmm_dir";
my $DIR_core_orthologs    = "$workdir/core_orthologs";
my $DIR_blast_db          = "$workdir/blast_dir";
my $DIR_HaMStR_tmp        = "$workdir/HaMStR_tmp";
my $DIR_HaMStR_res        = "$workdir/HaMStR_res";

system("mkdir -p $DIR_OrthoMCL_Fasta   ") unless -e $DIR_OrthoMCL_Fasta;
system("mkdir -p $DIR_OrthoMCL_DBmblast") unless -e $DIR_OrthoMCL_DBmblast;
system("mkdir -p $DIR_OrthoMCL_DBblast ") unless -e $DIR_OrthoMCL_DBblast;
system("mkdir -p $DIR_mclOutput        ") unless -e $DIR_mclOutput;
system("mkdir -p $DIR_group            ") unless -e $DIR_group;
system("mkdir -p $DIR_group_seq        ") unless -e $DIR_group_seq;
system("mkdir -p $DIR_group_aln        ") unless -e $DIR_group_aln;
system("mkdir -p $DIR_group_msf        ") unless -e $DIR_group_msf;
system("mkdir -p $DIR_hmmer            ") unless -e $DIR_hmmer;
system("mkdir -p $DIR_core_orthologs   ") unless -e $DIR_core_orthologs;
system("mkdir -p $DIR_blast_db         ") unless -e $DIR_blast_db;
system("mkdir -p $DIR_HaMStR_tmp       ") unless -e $DIR_HaMStR_tmp;
system("mkdir -p $DIR_HaMStR_res       ") unless -e $DIR_HaMStR_res;

my $softdir_OrthoMCL = "/mnt/ilustre/app/dna/orthologous_genes/orthomclSoftware-v2.0.9";
my $softdir_blast    = "/mnt/ilustre/app/dna/genome_annotation/annotation/blast+/bin";
my $softdir_mblastx  = "/mnt/ilustre/app/dna/genome_annotation/annotation/mblastx";
my $softdir_mcl_bin  = "/mnt/ilustre/app/rna/bin";
my $softdir_MAFFT    = "/mnt/ilustre/app/dna/multi_alignment/bin";
my $softdir_Gblocks  = "/mnt/ilustre/app/dna/multi_alignment/Gblocks_0.91b/Gblocks";
my $softdir_hmmer    = "/mnt/ilustre/app/rna/bin";
my $softdir_HaMStR   = "/mnt/ilustre/app/dna/orthologous_genes/hamstr.v13.1/bin";

my %Seqs;
my %Species_Retained;
my @Refer_Species;
my @Search_Species;
my $inflation_value_used = [split(/,/,$inflation_value_string)]->[0];
my $currentdir = getcwd;

&load_seqlist( $seqlist, \%Seqs ,\@Refer_Species, \@Search_Species, \%Species_Retained);

if($orthomcl){
	chdir $workdir;
	&step1_prepare_orthomcl_config($pjname);
	&step2_orthomclInstallSchema;
	&step3_orthomclAdjustFasta( \%{$Seqs{"orthomcl"}}, $DIR_OrthoMCL_Fasta );
	&step4_orthomclFilterFasta;
	&step5_All_to_All_BLAST;
	&step6_orthomclBlastParser;
	&step7_orthomclLoadBlast;
	&step8_orthomclPairs;
	&step9_orthomclDumpPairsFiles;
	chdir $currentdir;
}
&step10_mclCluster($inflation_value_string);
&step11_convert_mcl_to_groups($inflation_value_string);
&step12_obtain_1to1_OGs(\%{$Seqs{"orthomcl"}},$inflation_value_string,\$inflation_value_used);
&step13_abstract_seq_for_groups(\%{$Seqs{"orthomcl"}},$inflation_value_used);
&step14_multiple_alignment_for_each_group($inflation_value_used);
&step15_hmmbuild($inflation_value_used);
&step16_construct_blastdb_and_core_ortholog_set(\%{$Seqs{"orthomcl"}},$pjname,$inflation_value_used);
&step17_run_HaMStR(\%{$Seqs{"orthomcl"}},\%{$Seqs{"search"}},\@Refer_Species,$pjname,$inflation_value_used);
&step18_abstract_homologous_genes(\%{$Seqs{"search"}}, \@Refer_Species, \@Search_Species, $pjname, $inflation_value_used,$percents_string);

#####################################################
#####################################################
sub load_seqlist {
    my ( $seqlist, $hash, $refer ,$search, $retained_species) = @_;
    open FIN, "<" . $seqlist;
    while (<FIN>) {
        chomp;
        my @p = split /\t/, $_;
		next if $_ =~ /^#/;
        die "Error, seq list file must have three cols!\n"          if scalar @p < 2;
        die "species name should not be emplty!!\n"                 unless $p[0];
		die "type should not be empty!!\n"                          unless $p[1]; 
        die "file path should not be emplty!!\n"                    unless $p[2];
		die "type should be 'orthomcl' or 'search'!!\n"             if ($p[1] ne "orthomcl" && $p[1] ne "search");
        die "species name of type $p[1] should not be repeated!!\n" if exists $$hash{ $p[0] }{$p[1]};
        $p[2] = File::Spec->rel2abs( $p[2] );
        die "file of $p[0] does not exists!!\n" unless -e $p[2];

        $$hash{ $p[1] }{ $p[0] } = $p[2];
		push @{$refer},$p[0] if $p[1] eq "orthomcl";
		push @{$search},$p[0] if $p[1] eq "search";
    }
    close FIN;
	if($retainspeclist){
		open  FIN, "$retainspeclist";
		while(my $spec = <FIN>){
			chomp($spec);
			next unless exists $$hash{"search"}{$spec};
			$$retained_species{$spec}++;
		}
		close FIN;
	}
}

sub step1_prepare_orthomcl_config {
    my ($pjname) = @_;
    &print_time("###### step1_prepare_orthomcl_config");

    my $db_name  = "orthomcl_$pjname"; #数据库名
    my $db_host  = "20.1.1.220";    #主机名
    my $db_port  = '3306';         #端口号
	my $username = "stacks";         #用户名
    my $password = "stacks";        #密码
	my $dsn = "dbi:mysql:database=${db_name};hostname=${db_host};port=${db_port}";  #数据源

	my $drh = DBI->install_driver("mysql");
#	if (my $rc = $drh->func("dropdb",   $db_name, $db_host, $username, $password, "admin")){
#		print "drop database '", $db_name, "' successfully!\n";
#	}else{
#		print "database '", $db_name, "' does not exists!\n";
#	}

    #创建数据库$db_name
    if (my $rc = $drh->func("createdb", $db_name, $db_host, $username, $password, "admin")){
        print "create database '", $db_name, "' successfully!\n";
    }
    else {
        print "failed to create database '", $db_name, "'!\n";
    }

	###################
    my $config_file = "$workdir/orthomcl.config";
    copy("$softdir_OrthoMCL/doc/OrthoMCLEngine/Main/orthomcl.config.template",$config_file);
	`sed -i 's/dbConnectString=dbi:mysql:orthomcl:3307/dbConnectString=dbi:mysql:$db_name:20.1.1.220:3306/' $config_file`;
    `sed -i 's/dbLogin=/dbLogin=stacks/' $config_file`;
    `sed -i 's/dbPassword=/dbPassword=stacks/' $config_file`;
}

sub step2_orthomclInstallSchema {
	&print_time("###### step2_orthomclInstallSchema");
    my $config_file = "$workdir/orthomcl.config";
    my $config_log  = "$workdir/orthomcl.config.log";
    my $cmd = "$softdir_OrthoMCL/bin/orthomclInstallSchema $config_file $config_log";
	&print_time($cmd);
	system $cmd;
}

sub step3_orthomclAdjustFasta {
	&print_time("###### step3_orthomclAdjustFasta");
    my ( $hash, $outdir ) = @_;
    foreach my $species ( sort keys %{$hash} ) {
        my $seq_file = $$hash{$species};
		my $cmd_1 = "$softdir_OrthoMCL/bin/orthomclAdjustFasta $species $seq_file 1";
		&print_time($cmd_1);
		system $cmd_1;

		my $cmd_2 = "mv $workdir/$species.fasta $outdir";
		#my $cmd_2 = "mv $species.fasta $outdir";
		&print_time($cmd_2);
		system $cmd_2;
    }
}

sub step4_orthomclFilterFasta {
	&print_time("###### step5_orthomclFilterFasta");

    my $cmd = "$softdir_OrthoMCL/bin/orthomclFilterFasta $DIR_OrthoMCL_Fasta 10 20";
	&print_time($cmd);
	system $cmd;
}

sub step5_All_to_All_BLAST {
	&print_time("###### step5_All_to_All_BLAST");
	copy( "$softdir_mblastx/key.lic", "$workdir/key.lic" );
    copy( "$workdir/goodProteins.fasta", $DIR_OrthoMCL_DBmblast );
    my $cmd_1 = "$softdir_mblastx/mhashgen -d $DIR_OrthoMCL_DBmblast -f goodProteins -s goodProteins.fasta -T N";
    &print_time($cmd_1);
	system $cmd_1;

    my $cmd_2 = "$softdir_mblastx/mneighborgen -d $DIR_OrthoMCL_DBmblast";
    &print_time($cmd_2);
	system $cmd_2;

    my $cmd_3 = "$softdir_mblastx/mblastp -d $DIR_OrthoMCL_DBmblast -f S -c goodProteins -q $workdir/goodProteins.fasta -T 16 -M 112 -o $workdir/all_VS_all.blastout";
    &print_time($cmd_3);
	system $cmd_3;

	&print_time("convert mblast to m8 format");
	&mblast_to_m8_format("$workdir/all_VS_all.blastout");

	unlink("$workdir/key.lic");
}

sub mblast_to_m8_format {
    my ($file) = @_;
    my $file_m8 = "$file.m8";
    open FIN, "$file";
    open FOUT, ">$file_m8";
	my %hash;
	while (<FIN>) {
        chomp;
        next if (/^Query Id/);
        next if /^\s*$/;
        next if /^\s*Reference Database/;
        next if /^\s*Number of Sequences/;
        my @p = split /\t/, $_;
        print FOUT $p[0], "\t", $p[1], "\t", $p[4], "\t", $p[5], "\t", $p[6],
          "\t", $p[7], "\t", $p[8], "\t", $p[9], "\t", $p[10], "\t", $p[11],
          "\t", $p[2], "\t", $p[3], "\n";
    }
    close FIN;
    close FOUT;
	system("sort -k1,1 -k2,2 $file_m8 -o $file_m8");
}

sub step6_orthomclBlastParser {
	&print_time("###### step6_orthomclBlastParser");
    my $cmd = "$softdir_OrthoMCL/bin/orthomclBlastParser $workdir/all_VS_all.blastout.m8 $DIR_OrthoMCL_Fasta > $workdir/similarSequences.txt 2>/dev/null";
	&print_time($cmd);
	system($cmd);
}

sub step7_orthomclLoadBlast {
	&print_time("###### step7_orthomclLoadBlast");
    my $cmd = "$softdir_OrthoMCL/bin/orthomclLoadBlast $workdir/orthomcl.config $workdir/similarSequences.txt";
	&print_time($cmd);
	system $cmd;
}

sub step8_orthomclPairs{
	&print_time("###### step8_orthomclPairs");
	my $cmd = "$softdir_OrthoMCL/bin/orthomclPairs $workdir/orthomcl.config $workdir/pairs.log cleanup=yes";
	&print_time($cmd);
	system $cmd;
}

sub step9_orthomclDumpPairsFiles{
	&print_time("###### step9_orthomclDumpPairsFiles");
	my $cmd = "$softdir_OrthoMCL/bin/orthomclDumpPairsFiles $workdir/orthomcl.config";
	&print_time($cmd);
	system $cmd;


    my $db_name  = "orthomcl_$pjname"; #数据库名
    my $db_host  = "20.1.1.220";    #主机名
    my $db_port  = '3306';         #端口号
    my $username = "stacks";         #用户名
    my $password = "stacks";        #密码
	my $drh = DBI->install_driver("mysql");
	if (my $rc = $drh->func("dropdb",   $db_name, $db_host, $username, $password, "admin")){
		print "drop database '", $db_name, "' successfully!\n";
	}
}

sub step10_mclCluster{
	my ($inflation_value_string) = @_;
	&print_time("###### step10_mclCluster");
	
	my $abc_raw = "$workdir/mclInput";

	my $pm=new Parallel::ForkManager(5);
	my @inflation_values = split /,/,$inflation_value_string;
	foreach my $i(@inflation_values){
		$pm->start and next;
		my $cmd = "$softdir_mcl_bin/mcl $abc_raw --abc -I $i -o $DIR_mclOutput/$i\_mclout &> $DIR_mclOutput/$i\_log";
		&print_time($cmd);
		system($cmd);
		unlink("$DIR_mclOutput/$i\_log");
		$pm->finish;
	}
	$pm->wait_all_children;
}

sub step11_convert_mcl_to_groups{
	my ($inflation_value_string) = @_;
	&print_time("###### step11_convert_mcl_to_groups");

	my @inflation_values = split /,/,$inflation_value_string;
	foreach my $i(@inflation_values){
		my $cmd = "$softdir_OrthoMCL/bin/orthomclMclToGroups OG 1 < $DIR_mclOutput/$i\_mclout > $DIR_group/$i\_groups";
		&print_time($cmd);
		system($cmd);
	}
}

sub step12_obtain_1to1_OGs{
	my ($hash,$inflation_value_string,$inflation_value) = @_;
	&print_time("###### step12_obtain_1to1_OGs");
	
	my @inflation_values = split /,/,$inflation_value_string;
	my $max_num = 0;
	foreach my $i(@inflation_values){
		my $group_num = 0;
		open  FIN,  "$DIR_group/$i\_groups";
		open  FOUT,">$DIR_group/$i\_groups.1to1";
		while(my $line = <FIN>){
			chomp($line);
			my ($id,@elements) = split / /,$line;
			$id =~ s/://;
			my %num;
			foreach my $ele (@elements){
				my ($species,$gene) = split /\|/,$ele;
				$num{$species}++;
			}
			if(sum(values %num) == scalar(keys %num) && scalar(keys %num) == scalar(keys %{$hash})){
				print FOUT $line,"\n";
				$group_num++;
			}
		}
		close FOUT;
		close FIN;
		if($group_num > $max_num){
			$max_num = $group_num;
			$$inflation_value = $i;
		}
	}
}

####step13
sub step13_abstract_seq_for_groups{
	my ($hash,$inflation_value_string) = @_;
	my @inflation_values = split /,/,$inflation_value_string;
	&print_time("###### step13_abstract_seq_for_groups");

	my %Spec_Seqs;
	open FIN,"$workdir/goodProteins.fasta";
	$/ = "\n>";
	while(my $line = <FIN>){
		my ($header,@seqs) = split /\n/,$line;
		$header =~ s/>//g;
		my $sequence = join("",@seqs);
		$sequence =~ s/>//g;
		$Spec_Seqs{$header} = $sequence;
	}
    $/ = "\n";
    close FIN;

    my $pm = new Parallel::ForkManager(5);
    foreach my $i (@inflation_values) {
        $pm->start and next;
        open FIN, "$DIR_group/$i\_groups.1to1";
        while (<FIN>) {
			chomp;
            my ( $id, @p ) = split / /, $_;
            $id =~ s/://;
            system("mkdir -p $DIR_group_seq/$i") unless -e "$DIR_group_seq/$i";
            open FOUT, ">$DIR_group_seq/$i/$id.fasta";
            foreach my $gene (@p) {
                print FOUT ">$id|$gene\n", $Spec_Seqs{$gene}, "\n";
            }
            close FOUT;
        }
        close FIN;
        $pm->finish;
    }
    $pm->wait_all_children;
}

sub step14_multiple_alignment_for_each_group{
    my ( $inflation_value_string ) = @_;
    &print_time("###### step14_multiple_alignment_for_each_group");

    my $pm = new Parallel::ForkManager(8);
    my @inflation_values = split /,/, $inflation_value_string;
    foreach my $i (@inflation_values) {
		system("mkdir -p $DIR_group_aln/$i") unless -e "$DIR_group_aln/$i";
		foreach my $prot_file(glob("$DIR_group_seq/$i/*.fasta")){
			$pm->start and next;
			my $name = basename($prot_file);
			$name =~ s/.fasta$//;
			my $prot_file_mafft = "$DIR_group_aln/$i/$name.aln";
			system("$softdir_MAFFT/einsi --thread 3 --quiet $prot_file > $prot_file_mafft");
			$pm->finish;
		}
    }
    $pm->wait_all_children;
}

sub step15_hmmbuild{
	my ( $inflation_value_string ) = @_;
	&print_time("###### step15_hmmbuild");

	my $pm = new Parallel::ForkManager(10);
	my @inflation_values = split /,/, $inflation_value_string;
	foreach my $i (@inflation_values) {
		system("mkdir -p $DIR_hmmer/$i") unless -e "$DIR_hmmer/$i";
		system("mkdir -p $DIR_group_msf/$i") unless -e "$DIR_group_msf/$i";
		foreach my $msa_file(glob("$DIR_group_aln/$i/*.aln")){
			$pm->start and next;
			my $name = basename($msa_file);
			$name =~ s/.aln$//;
			my $file_hmm = "$DIR_hmmer/$i/$name.hmm";
			my $file_log = "$DIR_hmmer/$i/$name.log";
			my $file_msf = "$DIR_group_msf/$i/$name.msf";
			my $cmd = "$softdir_hmmer/hmmbuild -n $name -o $file_log -O $file_msf $file_hmm $msa_file";
			system($cmd);
			unlink($file_log);
			$pm->finish;
		}
	}
    $pm->wait_all_children;
}

sub step16_construct_blastdb_and_core_ortholog_set{
	my ( $hash, $pjname, $inflation_value_string ) = @_;
	&print_time("###### step16_construct_blastdb_and_core_ortholog_set");

	foreach my $species(sort keys %{$hash}){
		my $dir  = "$DIR_blast_db/$species";
		my $file = "$DIR_blast_db/$species/$species\_prot.fa";
		my $name = "$DIR_blast_db/$species/$species\_prot";
		mkdir $dir unless -e $dir;
		copy($$hash{$species},$file);
		system("$softdir_blast/makeblastdb -in $file -input_type fasta -dbtype prot -parse_seqids -out $name &> /dev/null");
		#system("/mnt/lustre/share/apps/blast-2.2.26/bin/formatdb -i $file -n $name -t $name -o T -l $name.log -p T");
	}

	my @inflation_values = split /,/, $inflation_value_string;
	foreach my $i (@inflation_values) {
		my $DIR_core_orthologs_set = "$DIR_core_orthologs/$i/$pjname";
		system("mkdir -p $DIR_core_orthologs_set") unless -e $DIR_core_orthologs_set;
		system("mkdir -p $DIR_core_orthologs_set/hmm_dir");
		system("mkdir -p $DIR_core_orthologs_set/aln_dir");
		system("mkdir -p $DIR_core_orthologs_set/aln_dir/fa_dir");
		system("mkdir -p $DIR_core_orthologs_set/aln_dir/msf_dir");

		system("cat $DIR_group_seq/$i/*.fasta > $DIR_core_orthologs_set/$pjname.fa");
		system("cp  $DIR_hmmer/$i/*.hmm         $DIR_core_orthologs_set/hmm_dir/");
		system("cp  $DIR_group_aln/$i/*.aln     $DIR_core_orthologs_set/aln_dir/fa_dir");
		system("cp  $DIR_group_msf/$i/*.msf     $DIR_core_orthologs_set/aln_dir/msf_dir");
	}
}

sub step17_run_HaMStR{
	my ( $hash_1, $hash, $refer, $pjname, $inflation_value_string ) = @_;
	&print_time("###### step17_run_HaMStR");

	my @inflation_values = split /,/, $inflation_value_string;
	my $refer_strings = join(",", @{$refer});

	my $pm = new Parallel::ForkManager(5);
	foreach my $i (@inflation_values) {
		my $DIR_core_orthologs_i = "$DIR_core_orthologs/$i";
		foreach my $species(keys %$hash){
			$pm->start and next;
			my $file = $$hash{$species};
			my $outdir_tmp = "$DIR_HaMStR_tmp/$i/$species";
			my $outdir_res = "$DIR_HaMStR_res/$i/";
			system("mkdir -p $outdir_tmp") unless -e $outdir_tmp;
			system("mkdir -p $outdir_res") unless -e $outdir_res;
			#if(@{$refer} ~~ /$species/){
			if(exists $$hash_1{$species}){
				my $refer_prime = $$refer[0];
				foreach my $OG_file(glob("$DIR_group_seq/$i/*.fasta")){
					my $name = basename($OG_file);
					my $dir = "$DIR_HaMStR_tmp/$i/$species/fa_dir_$species\_$pjname\_$refer_prime\_strict";
					mkdir $dir unless -e $dir;
					my $SS;
					my $HH;

					open FIN , "$OG_file";
					open FOUT, ">$dir/$name";
					$/="\n>";
					while(my $line = <FIN>){
						chomp($line);
						my ($header,@seqs) = split /\n/,$line;
						$header =~ s/>//;
						my $sequence = join("",@seqs);
					print FOUT ">$header\n$sequence\n";

						my ($og,$spec,$gene) = split /\|/,$header;
						if($spec eq $species){
							$HH = "$og|$spec|$spec|$gene|1";
							$SS = $sequence;
						}
					}
					$/="\n";
					print FOUT ">$HH\n$SS\n";
					close FIN;
					close FOUT;
				}
			}else{
				my $cmd = "$softdir_HaMStR/hamstr_modified.pl -sequence_file=$file -taxon=$species -hmmset=$pjname -hmmpath=$DIR_core_orthologs_i -refspec=$refer_strings -blastpath=$DIR_blast_db -strict -representative -rbh -eval_hmmer=0.00001 -outpath=$outdir_tmp -cpu 5 -force &> $outdir_tmp/log.txt";
				&print_time($cmd);
				system("$cmd");
			}
			$pm->finish;
		}
	}
	$pm->wait_all_children;
}

sub step18_abstract_homologous_genes{
	my ( $hash, $refer, $search, $pjname, $inflation_value_string, $percents_string ) = @_;
	&print_time("###### step18_abstract_homologous_genes");

	system("rm -rf $DIR_HaMStR_res") if -e $DIR_HaMStR_res;

	my @inflation_values = split /,/, $inflation_value_string;

	my %H;
	my $refer_prime = $$refer[0];
	foreach my $i(@inflation_values){
		foreach my $species(sort keys %{$hash}){
			my @Files = glob("$DIR_HaMStR_tmp/$i/$species/fa_dir_$species\_$pjname\_$refer_prime\_strict/*");
			foreach my $file(@Files){
				my $OG_Name = basename($file);
				$OG_Name =~ s/\.fasta//;
				open  FIN, "<$file";
				$/="\n>";
				while(my $line = <FIN>){
					chomp($line);
					my ($header,@seqs) = split /\n/,$line;
					$header =~ s/>//;
					my $sequence = join("",@seqs);
					my @p = split /\|/,$header;
					if(scalar @p == 3){
						my ($OG,$species,$gene) = @p;
						$H{$i}{$OG}{$species} = ">$species\__$gene\n$sequence\n";
					}else{
						my ($OG,$ref_spec,$species,$gene,$num) = @p;
						$H{$i}{$OG}{$species} = ">$species\__$gene\n$sequence\n";
					}
				}
				$/="\n";
				close FIN;
			}
		}
	}

	foreach my $i(sort keys %H){
		my @Percents = split /,/,$percents_string;
		foreach my $OG(sort keys %{$H{$i}}){
			my $flag=0;
			foreach my $spec(keys %Species_Retained){
				unless(exists($H{$i}{$OG}{$spec})){
					$flag = 0;
					last;
				}else{
					$flag = 1;
				}
			}
			if($flag){
				system "mkdir -p $DIR_HaMStR_res/$i/retained" unless -e "$DIR_HaMStR_res/$i/retained";
				open  FOUT,">$DIR_HaMStR_res/$i/retained/$OG.fa";
				foreach my $species(@{$search}){
					if(exists $H{$i}{$OG}{$species}){
						print FOUT $H{$i}{$OG}{$species};
					}
				}
				close FOUT;
			}
				
			foreach my $perc(@Percents){
				my $num=0;
				foreach my $spec(keys %{$H{$i}{$OG}}){
					$num++ if exists $$hash{$spec};
				}
				my $total = scalar @{$search};
				my $ratio = $num/$total;
				#print $OG,"\t",$num,"\t",$total,"\t",$ratio,"\t",$perc,"\n";

				if($ratio >= $perc){
					system "mkdir -p $DIR_HaMStR_res/$i/$perc" unless -e "$DIR_HaMStR_res/$i/$perc";
					open  FOUT,">$DIR_HaMStR_res/$i/$perc/$OG.fa";
					foreach my $species(@{$search}){
						if(exists $H{$i}{$OG}{$species}){
							print FOUT $H{$i}{$OG}{$species};
						}
					}
					close FOUT;
				}
			}
		}
		######################
		my @elements;
		if(-e "$DIR_HaMStR_res/$i/retained"){
			@elements = (@Percents,"retained");
		}else{
			@elements = @Percents;
		}

		foreach my $ele(@elements){
			my %HHH;
			my %SPEC;
			foreach my $file(glob("$DIR_HaMStR_res/$i/$ele/*.fa")){
				my $OG = basename($file);
				$OG    =~ s/.fa$//;
				open  FIN,$file;
				while(my $line = <FIN>){
					chomp($line);
					next unless $line =~ /^>/;
					$line =~ s/^>//;
					my ($species,$gene) = split /__/,$line;
					$SPEC{$species}++;
					$HHH{$OG}{$species} = $gene;
				}
				close FIN;
			}
			open  FOUT,">$DIR_HaMStR_res/$i/$ele/Matrix.txt";
			print FOUT "OGs","\t",join("\t",sort keys %SPEC),"\n";
			foreach my $OG(sort keys %HHH){
				print FOUT $OG;
				foreach my $species(sort keys %SPEC){
					if(exists($HHH{$OG}{$species})){
						print FOUT "\t",$HHH{$OG}{$species};
					}else{
						print FOUT "\t";
					}
				}
				print FOUT "\n";
			}
			close FOUT;
		}
	}
}

sub sum {
    my @elements = @_;
    my $sum      = 0;
    foreach my $ele (@elements) {
        $sum += $ele;
    }
    return $sum;
}

sub getTime() {
    (
        my $sec, my $min, my $hour, my $day, my $mon, my $year, my $weekday,
        my $yeardate, my $savinglightday
    ) = ( localtime(time) );
    $sec  = ( $sec < 10 )  ? "0$sec"  : $sec;
    $min  = ( $min < 10 )  ? "0$min"  : $min;
    $hour = ( $hour < 10 ) ? "0$hour" : $hour;
    $day  = ( $day < 10 )  ? "0$day"  : $day;
    $mon = ( $mon < 9 ) ? "0" . ( $mon + 1 ) : ( $mon + 1 );
    $year += 1900;
    my $now = "$year-$mon-$day $hour:$min:$sec";
    return $now;
}

sub print_time {
    my @lines = @_;
    foreach my $line (@lines) {
        print "[", &getTime . "] : " . $line . "\n";
    }
}
