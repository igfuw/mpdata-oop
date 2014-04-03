#

length($0) != 0 && ($2 != "started" ) { 
  fc = substr($1,1,1); 
  if (fc == "[" || fc == "#") 
  { 
    ll=" "$0
  } 
  else 
  {
    printf("%s\n%s", ll, $0); 
    ll="";
  } 
} 

END {
  print
}
