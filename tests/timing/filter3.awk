# hides load-only records and appends the loadtime column to calc runs

BEGIN {FS=";"}                                             
{                                                          
 if (NR % 2 == 1) offset=$13;                             
 else print $0 "loadtime=;" offset;                      
}                                                          
