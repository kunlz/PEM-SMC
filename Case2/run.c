#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <dirent.h>
#include <fcntl.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>
#include <sys/types.h> 
#include <math.h>
 


int main(int argc, char **argv)
{
    
//  run the jobclm.single.csh file
    system(" chmod  u+x /group_homes/lzu_public/home/u120220909911/Summer/Last/LE/*.csh");
    system("/group_homes/lzu_public/home/u120220909911/Summer/Last/LE/jobclm.csh");
}

