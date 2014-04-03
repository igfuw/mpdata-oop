#  

{                               
  if (index($0,"#") != 0)       
  {                              
    l=l $2 "=;" $4 "%;";         
  } else if (NR==1) {
    printf("%s", $0)
  } else {                       
    printf("%s", l "\n" $0);     
    l="";                        
  }                              
} 
