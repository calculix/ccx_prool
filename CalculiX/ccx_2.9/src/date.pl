#!/usr/bin/perl

chomp($date=`date`);

# inserting the date into ccx_2.9.c

@ARGV="ccx_2.9.c";
$^I=".old";
while(<>){
    s/You are using an executable made on.*/You are using an executable made on $date\\n");/g;
    print;
}

system "rm -f ccx_2.9.c.old";
