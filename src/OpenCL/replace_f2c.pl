#!/usr/bin/perl

while (<>) {
while (/(cov|diu)(\d)\(i,j,k\)/) {
my $v=$1;
my $np1=($2)-1;

s/(cov|diu)(\d)\(i,j,k\)/${1}_ijk.s$np1/;
#print $_;
};

/(cov|diu)(\d)\(i\+1,j,k\)/ && do {
my $v=$1;
my $np1=($2)-1;

s/(cov|diu)(\d)\(i\+1,j,k\)/${1}_ijk_p1.s$np1/;
#print $_;
};

/(cov|diu)(\d)\(i,j\+1,k\)/ && do {
my $v=$1;
my $np1=($2)-1;

s/(cov|diu)(\d)\(i,j\+1,k\)/${1}_ijk_p1.s$np1/;
#print $_;
};

/(cov|diu)(\d)\(i,j,k\+1\)/ && do {
my $v=$1;
my $np1=($2)-1;

s/(cov|diu)(\d)\(i,j,k\+1\)/${1}_ijk_p1.s$np1/;
#print $_;
};

s/dx1\((.+?)\)/dx1[$1+1]/g;
s/dy1\((.+?)\)/dx1[$1]/g;
s/dzs\((.+?)\)/dx1[$1+1]/g;
s/dzn\((.+?)\)/dx1[$1+1]/g;

s/\s\&\s*$//;
s/([\.\)\]\d])\s*$/$1;\n/;
print $_;
}