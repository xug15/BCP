open I, "<$ARGV[0]";

@otuname=qw//;
while(<I>){
#print "$_";
chomp;
if($_=~/>(.*?)$/){
    $name=$1;

push(@otuname, $name);

}

}


close I;
$n=0;
$total=length(@otuname);
print "lenth\t$#otuname\n";

foreach(@otuname){
    print "$n\n";
print "$otuname[$n]\n";

#push (@newlist, 1, $n);
$info=join "\t", @newlist;
print "$info\n";
$n++;
}
push @newarray, 'TT' foreach (1..$total);
print "@newlist\n";

@newlist = (@genotypes, (0) x 20);
print "@newlist\n";
