#!/usr/bin/perl 
use strict;
use 5.012;

use File::Copy;
use F95Normaliser qw( normalise_F95_src );

my $f95src_path='../RefactoredSources';
my $csrc_path='../RefactoredSources/PostCPP';
my $kernels_path='./Kernels';
my $stubs_path='./KernelStubs';
my $wrapper_templ_path='./WrapperTemplates';

@ARGV || die "Please specify the kernel name.\n";
my $kname=$ARGV[0];
$kname=~s/\..+$//;


say "Copying C code to ${kname}_kernel.cl";
copy "$csrc_path/$kname.c","$kernels_path/${kname}_kernel.cl";

open my $F95SRC, '<', "$f95src_path/$kname.f95" || die "$!";

# First create the template for the OpenCL module
open my $F95MODTMPL, '>', "$wrapper_templ_path/module_${kname}_ocl_TEMPL_V2.f95" || die "$!";
my @orig_lines=<$F95SRC>;
close $F95SRC;
my $src_lines = normalise_F95_src(\@orig_lines);

#use Data::Dumper;
#die Dumper($src_lines);

my $decl=2;
my $sub_args='';
for my $line (@{$src_lines} ) {
# cheap check to insert kernel call
# Detect variable declarations
    if ($line!~/^\s*!/ && $line !~/^\s*$/ && $line=~/::/) {$decl=1};
# Detect the line after the last declaration    
    if ($line!~/^\s*!/ && $line !~/^\s*$/ && $line!~/::/ && $decl==1) {
        $decl=0;
    }
# Detect the end of the subroutine    
    if ($line=~/end\s+subroutine/ && $decl==3) {
        $decl=2;
    }
# Detect the beginning of the module, and transform into OpenCL module    
    $line=~/^\s*module\s+(\w+)/ && do { 
        my $modname=$1;
        print $F95MODTMPL "module ${modname}_ocl\n";
        print $F95MODTMPL "    use module_LES_conversions\n";
        print $F95MODTMPL "    use $modname\n";

        next;
    };
# Detect the end of the module, and transform into OpenCL module    
    $line=~/^\s*end\s+module\s+(\w+)/ && do {
        my $modname = $1;
        my $tline=$line;
        $tline=~s/$modname/${modname}_ocl/;
        if ($decl==3) {
        warn "Did not find 'end subroutine', inserting ad-hoc!\n";
        say $F95MODTMPL '    end subroutine';
        }
        say $F95MODTMPL $tline;
        next;
    };
# Detect subroutine and transform into OpenCL subroutine
    $line=~/^\s*subroutine\s+(\w+)\s*\((.+)\)/ && do {
        my $subname = $1;
        $sub_args=$2;
        my $tline=$line;
        $tline=~s/$subname/${subname}_ocl/;
        say $F95MODTMPL $tline;
        next;
    };
# Ditto at end
    $line=~/end\s+subroutine\s+(\w+)/ && do {
        my $subname = $1;
        my $tline=$line;
        $tline=~s/$subname/${subname}_ocl/;
        say $F95MODTMPL $tline;
        next;
    };
# Insert pragmas and subroutine call after declarations    
    if ($decl==0) {
        print $F95MODTMPL "    !\$ACC Kernel(${kname}_kernel), GlobalRange(1), LocalRange(1) \n";
        print $F95MODTMPL "    call $kname($sub_args)\n";
        print $F95MODTMPL "    !\$ACC End Kernel\n";
        $decl=3;
        next;
    } elsif ($decl<3) {
    say $F95MODTMPL $line;
    }

}
close $F95MODTMPL;

# Next create a Fortran stub for the kernel, code very similar to the above
open my $F95KSTUB, '>', "$stubs_path/${kname}_kernel.f95" || die "$!";
$decl=2;
for my $line (@{$src_lines} ) {
# cheap check to insert kernel call
    if ($line!~/^\s*!/ && $line !~/^\s*$/ && $line=~/::/) {$decl=1};
    if ($line!~/^\s*!/ && $line !~/^\s*$/ && $line!~/::/ && $decl==1) {
        $decl=0;
    } 
    if ($line=~/end\s+subroutine/ && $decl==3) {
        $decl=2;
    }
    
    $line=~/^\s*module\s+(\w+)/ && do { 
        my $modname=$1;

        print $F95KSTUB "module ${modname}_kernel\n";
        next;
    };
    $line=~/^\s*end\s+module\s+(\w+)/ && do {
        my $modname = $1;
        my $kline=$line;
        $kline=~s/$modname/${modname}_kernel/;
        if ($decl==3) {
        warn "Did not find 'end subroutine', inserting ad-hoc!\n";
        say $F95KSTUB '    end subroutine';
        }
        say $F95KSTUB $kline;
        next;
    };

    $line=~/subroutine\s+(\w+)/ && do {
        my $subname = $1;
        my $kline=$line;
        $kline=~s/$subname/${subname}_kernel/;
        say $F95KSTUB $kline;
        next;
    };
    if ($decl==0) {
        $decl=3;
        next;
    } elsif ($decl<3) {
        say $F95KSTUB $line;
    }

}
close $F95KSTUB;

=info
This script does three things:

1. Copy the C code generated with F2C_ACC into Kernels and rename to .cl. Later, it will also clean up the C code so that it's valid OpenCL.
2. Create a F95 OpenCL wrapper template and put it in WrapperTemplates
3. Create a F95 kernel stub an put it in KernelStubs

This is a very cheap, non-robust way to generate the OpenCL Fortran wrapper template and the Fortran kernel stub.
The script is dumb, it will do the transformation for every routine in the original module.  
=cut
