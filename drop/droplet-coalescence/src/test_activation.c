#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "persistent_activation.h"
int main(void){
  int first[4]={-1,-1,-1,-1};double h[4]={.030,.031,.02,.04};
  assert(activation_update(4,first,h,.03,0)==2);
  assert(first[0]==-1 && first[1]==0 && first[2]==-1 && first[3]==0);
  int ids[3]={1,3,2};double coeff[3]={-1.,.7,.3};
  assert(!activation_face_allowed(3,ids,coeff,first));
  h[1]=.005;h[2]=.031;
  assert(activation_update(4,first,h,.03,1)==1);
  assert(first[1]==0 && first[2]==1);
  assert(activation_face_allowed(3,ids,coeff,first));
  h[2]=.001;assert(activation_update(4,first,h,.03,2)==0 && first[2]==1);
  double source[4]={0};for(int k=0;k<3;k++)source[ids[k]]-=coeff[k]*2.;
  assert(source[0]==0. && fabs(source[0]+source[1]+source[2]+source[3])<1e-14);
  coeff[2]=0.;first[2]=-1;assert(activation_face_allowed(3,ids,coeff,first));
  puts("threshold strictness, persistence, all-stencil masking, zero source and conservation PASS");return 0;
}
