#!/usr/bin/env perl

chomp($date=`date`);

# inserting the date into ccx_2.19.c

@ARGV="ccx_2.19.c";
$^I=".old";
while(<>){
    s/You are using an executable made on.*/You are using an executable made on $date\\n");/g;
    print;
}

# inserting the date into ccx_2.19step.c

@ARGV="ccx_2.19step.c";
$^I=".old";
while(<>){
    s/You are using an executable made on.*/You are using an executable made on $date\\n");/g;
    print;
}

# inserting the date into frd.c

@ARGV="frd.c";
$^I=".old";
while(<>){
    s/COMPILETIME.*/COMPILETIME       $date                    \\n\",p1);/g;
    print;
}

system "rm -f ccx_2.19.c.old";
system "rm -f ccx_2.19step.c.old";
system "rm -f frd.c.old";
