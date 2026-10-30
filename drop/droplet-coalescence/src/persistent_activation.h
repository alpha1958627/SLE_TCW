#ifndef PERSISTENT_ACTIVATION_H
#define PERSISTENT_ACTIVATION_H

/* A cell activates strictly above the threshold and never deactivates.
   With reconstructed AMR gradients, every nonzero stencil participant must
   be active before a face can inject noise into the adjoint divergence. */
static int activation_update(int n,int *first_step,const double *height,double threshold,int step){
  int added=0;
  for(int i=0;i<n;i++)if(first_step[i]<0 && height[i]>threshold){first_step[i]=step;added++;}
  return added;
}
static int activation_face_allowed(int ng,const int *ids,const double *g,const int *first_step){
  for(int k=0;k<ng;k++)if(g[k]!=0. && first_step[ids[k]]<0)return 0;
  return 1;
}
#endif
