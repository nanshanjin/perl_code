use strict;
open IN,"/home/Account/cuilefang/prjt/1118/group_info_n.txt";
my @infos;
while(<IN>){
    $_=~s/[\n\r]//g;
    next if /^PATIENT_ID/;
    my ($id,undef)=split/\t/,$_;
    push(@infos,$id);
    }
close IN;

open IN,"/home/Account/cuilefang/prjt/1118/RSEM_normd.txt";
my @nums;
my %hash;
while(<IN>){
    $_=~s/[\n\r]//g;
    if(/^Hugo_Symbol/){
        my ($Hugo_Symbol,$Entrez_Gene_Id,@ids)=split/\t/,$_;
        for(my $i=0;$i<=$#ids;$i++){
            $hash{$ids[$i]}=$i;
            }
        foreach my $keys (keys %hash){
            for(my $i=0;$i<=$#infos;$i++){
                if($keys eq $infos[$i]){
                    push(@nums,$hash{$keys});
                    }
                }
            }
        print join("\t","Hugo_Symbol","Entrez_Gene_Id",@infos),"\n";
        }else{
            my ($Hugo_Symbol,$Entrez_Gene_Id,@ids)=split/\t/,$_;
            my @tmps;
            for(my $i=0;$i<=$#nums;$i++){
                push(@tmps,$ids[$nums[$i]]);
                }
            print scalar @nums;
            print join("\t",$Hugo_Symbol,$Entrez_Gene_Id,@tmps),"\n";
        }
    }
close IN;
