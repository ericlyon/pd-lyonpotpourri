if (! -e "lyonpotpourri"){
    mkdir("lyonpotpourri");
}
while(<*>){
    chomp;
    if(/darwin$/ | /liblyonpotpourri.dylib/ | /help.pd/){
	`mv $_ lyonpotpourri`;
    }
}
`cp lpp-icon.pd lyonpotpourri`;
`cp lpp-meters.pd lyonpotpourri`;
`mv sound lyonpotpourri`;
