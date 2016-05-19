#include <stdio.h>

#define P0 0
#define P1 1.0
#define P2 19.0
//#define P3 74.0
#define P3 109.0
#define P4 129.0
//#define P5 184.0
#define P5 149.0
//#define P6 239.0
#define P6 215.0

#define BLDV 50
#define BLTHICK 10
#define DV .75
int main(void)
{

  int i,j;
  unsigned char pallete[768];
  long int palout;
  FILE*in;
  FILE*out;
  float temp;

  if( (in=fopen("john.pal","rb"))==NULL){
    printf("error: in\n");
    return(1);
  }
  if( (out=fopen("jon2.pal","wb"))==NULL){
    printf("error: out\n");
    return(1);
  }
  if( (out=fopen("jon2.rpal","wb"))==NULL){
    printf("error: rout\n");
    return(1);
  }

  for(i=0;i<768;i++){
    fread(&pallete[i],sizeof(unsigned char),1,in);
  }
 
  printf("%5s %3s %3s %3s\n","value","R","G","B");
  for(i=0;i<256;i++){
    printf("%5d %3d %3d %3d\n",i,(short)pallete[i*1],(short)pallete[i+256],(short)pallete[i+256*2]);
  }

  printf("now new-----\n");
  
  for(i=0;i<768;i++){
    pallete[i]=0;
  }

  // red
  for(i=0;i<256;i++){
    if(i<P4){
      pallete[i]=0;
    }
    if( (i>=P4)&&(i<=P5)){
      pallete[i]=(short)(((float)i-P4)*(255.0/(P5-P4)));
    }
    if(i>P5){
      pallete[i]=255;
    }
  }
  // green
  for(i=0;i<256;i++){
    if(i<P2){
      pallete[i+256]=0;
    }
    if( (i>=P2)&&(i<=P3)){
      pallete[i+256]=(short)(((float)i-P2)*(255.0/(P3-P2)));
    }
    if( (i>P3)&&(i<P5)){
      pallete[i+256]=255;
    }
    if( (i>=P5)&&(i<=P6)){
      pallete[i+256]=(short)(((float)i-P6)*(-255.0/(P6-P5)));
    }
    if( (i>=P6)&&(i<=255)){
      pallete[i+256]=(short)(((float)i-P6)*(255.0/(255.-P6)));
    }
  }
  // blue
  for(i=0;i<256;i++){
    if(i==P0) pallete[i+256*2]=0;
    if( (i>=1)&&(i<P2)){
      pallete[i+256*2]=(short)(((float)i-1.0)*((255.0-128.)/(P2-1.)))+128;
    }
    if(i==P2) pallete[i+256*2]=255; // above for some reason doesn't work here
    if( (i>P2)&&(i<P3)){
      pallete[i+256*2]=255;
    }
    if( (i>=P3)&&(i<=P4)){
      pallete[i+256*2]=(short)(((float)i-P4)*(-255.0/(P4-P3)));
    }
    if( (i>=P4)&&(i<=P6)){
      pallete[i+256*2]=0;
    }
    if( (i>=P6)&&(i<=255)){
      pallete[i+256*2]=(short)(((float)i-P6)*(255.0/(255.-P6)));
    }
  }
  /*
  for(i=0;i<256;i++){
    if(!(i%(BLDV))){
       for(j=0;j<BLTHICK;j++){
	 if(i+j+256*2<256*3){
	   if(DV>0){
	     pallete[i+j]=(short)((float)(pallete[i+j])*DV);
	   }
	   else pallete[i+j]=0;
	   if(DV>0){
	     pallete[i+j+256]=(short)((float)(pallete[i+j+256])*DV);
	   }
	   else pallete[i+j+256]=0;
	   if(DV>0){
	     pallete[i+j+256*2]=(short)((float)(pallete[i+j+256*2])*DV);
	   }
	   else pallete[i+j+256*2]=0;
	 }
       }
    }
  }
  */
  /*
  for(i=(P3+P2)/2-BLTHICK/2;i<(P3+P2)/2+BLTHICK/2;i++){
    pallete[i]/=DV;
    pallete[i+256]/=DV;
    pallete[i+256*2]/=DV;
  }
  
  for(i=(P4+P3)/2-BLTHICK/2;i<(P4+P3)/2+BLTHICK/2;i++){
    pallete[i]/=DV;
    pallete[i+256]/=DV;
    pallete[i+256*2]/=DV;
  }
  
  for(i=(P5+P4)/2-BLTHICK/2;i<(P5+P4)/2+BLTHICK/2;i++){
    pallete[i]/=DV;
    pallete[i+256]/=DV;
    pallete[i+256*2]/=DV;
  }
  for(i=(P6+P5)/2-BLTHICK/2;i<(P6+P5)/2+BLTHICK/2;i++){
    pallete[i]/=DV;
    pallete[i+256]/=DV;
    pallete[i+256*2]/=DV;
  }
  */
  for(i=(255+P6)/2-BLTHICK/2;i<(255+P6)/2+BLTHICK/2;i++){
    pallete[i]/=DV;
    pallete[i+256]/=DV;
    pallete[i+256*2]/=DV;
  }
  
  for(i=0;i<768;i++){
    fwrite(&pallete[i],sizeof(unsigned char),1,out);
  }

  printf("%5s %3s %3s %3s\n","value","R","G","B");
  for(i=0;i<256;i++){
    printf("%5d %3d %3d %3d\n",i,(short)pallete[i*1],(short)pallete[i+256],(short)pallete[i+256*2]);
  }
  
  
  return(0);
}
