#ifndef AMR_HEIGHT_H
#define AMR_HEIGHT_H
static double amr_height_lower=1e-6;
static long amr_transfer_limited=0;
/* Limited slopes, one common factor for all children. Child mean equals
   parent exactly; no volume is added and no child is height-clipped. */
static void amr_refine_positive (Point point,scalar s) {
  double base=s[],budget=fmax(base-amr_height_lower,0.),extent=0.;coord g;
  foreach_dimension(){g.x=(s[1]-s[-1])/8.;extent+=fabs(g.x);}
  double factor=extent>0.?fmin(1.,budget/extent):1.;
  if(factor<1.-1e-12&&is_local(cell))amr_transfer_limited++;
  foreach_child(){s[]=base;foreach_dimension()s[]+=factor*child.x*g.x;}
}
static double amr_cos4_point (double xx,double yy) {
  double r1=hypot(xx-(L0/2.-4.75),yy-L0/2.);
  double r2=hypot(xx-(L0/2.+4.75),yy-L0/2.);
  return .01+(r1<5.?pow(cos(pi*r1/10.),4):0.)+(r2<5.?pow(cos(pi*r2/10.),4):0.);
}
/* Exact restriction of the SAME finest-grid discrete initial condition,
   without ever allocating a full finest tree for the AMR simulation. */
static double amr_initial_average (double xx,double yy,double d,int maximum) {
  double rx1=fmax(fabs(xx-(L0/2.-4.75))-.5*d,0.);
  double rx2=fmax(fabs(xx-(L0/2.+4.75))-.5*d,0.);
  double ry=fmax(fabs(yy-L0/2.)-.5*d,0.);
  if(rx1*rx1+ry*ry>=25.&&rx2*rx2+ry*ry>=25.)return .01;
  double fine=L0/(1<<maximum);int m=(int)llround(d/fine);double sum=0.;
  for(int j=0;j<m;j++)for(int i=0;i<m;i++)sum+=amr_cos4_point(xx-.5*d+(i+.5)*fine,yy-.5*d+(j+.5)*fine);
  return sum/(m*(double)m);
}
static void amr_initialise (scalar h,scalar pressure,int minimum,int maximum,
                            double tolerance) {
  h.refine=amr_refine_positive;h.prolongation=refine_linear;
  h.restriction=restriction_volume_average;
  foreach(){h[]=amr_initial_average(x,y,Delta,maximum);pressure[]=0.;}
  boundary({h,pressure});
  scalar mark[];
  for(int target=minimum+1;target<=maximum;target++) {
    foreach() {
      coord g;foreach_dimension()g.x=(h[1]-h[-1])/8.;
      double worst=0.;
      for(int jy=-1;jy<=1;jy+=2)for(int ix=-1;ix<=1;ix+=2) {
        double actual=amr_initial_average(x+ix*Delta/4.,y+jy*Delta/4.,Delta/2.,maximum);
        worst=fmax(worst,fabs(actual-h[]-ix*g.x-jy*g.y));
      }
      mark[]=worst;
    }
    boundary({mark});refine(mark[]>tolerance&&level<target);amr_grid_sync();
    foreach(){h[]=amr_initial_average(x,y,Delta,maximum);pressure[]=0.;}
    boundary({h,pressure});
  }
}
#endif
