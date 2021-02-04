if (! -e "lyonpotpourri"){
    mkdir("lyonpotpourri");
}
while(<*>){
    chomp;
    if(/darwin$/ | /liblyonpotpourri.dylib/){
	`mv $_ lyonpotpourri`;
    }
}
