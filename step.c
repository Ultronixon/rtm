void step_forward(float * tgt_p, 
    float * src_p,
    float * v,
    int nz, int nx,
    float rz, float rx, float s)
/*< step forward >*/
{

  int iz;
  int ix;
  int ioff;
  float two=2.0;

  for (ix=1;ix<nx-1;ix++) {
    for (iz=1;iz<nz-1;iz++) {
      ioff=iz+ix*nz;
      tgt_p[ioff]=two*src_p[ioff]-tgt_p[ioff] +v[ioff]*(rz*(src_p[ioff+1]+src_p[ioff-1]) +
                  rx*(src_p[ioff+nz]+src_p[ioff-nz]) - s*src_p[ioff]);
    }
  }
}


