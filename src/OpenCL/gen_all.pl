my @templs = glob('WrapperTemplates/*');
for my $templ(@templs) {
system("./oclgen -o $templ");
}
system('./oclgen module_LES_combined_ocl_TEMPL_V2.f95 > Wrappers/module_LES_combined_ocl.f95');
