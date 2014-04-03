# choses the shortes run from each 3 runs

BEGIN { 
  FS=";"; 
  huge=1e6; 
  besttime=huge;
} 

{ 
   lines[NR%3]=$0; 
   if (besttime > $13) 
   { 
     besttime = $13; 
     bestid=NR%3; 
   } 
   if (NR%3 == 0) 
   { 
     print lines[bestid]; 
     besttime=huge; 
   } 
}                               
