#include "KFKSDS-steady-v2.h"

void KFKSDS_steady_v2eps (int *dim, double *sy, double *sZ, double *sZtZ, 
  double *sT, double *sH, double *sR, 
  double *V, double *sQ, double *sa0, double *sP0, 
  double *tol, int *maxiter, double *ksconvfactor,
  double *res)
{
  int i, ip1, n = dim[0], m = dim[2], //ir = dim[3], 
    convref, nmconvref, nm1 = n-1;
  
  //double v[n], f[n], invf[n], vof[n];
  std::vector<double> v(n), f(n), invf(n), vof(n);
  double sHsq = pow(*sH, 2), msHsq = -1.0 * sHsq,
    epshat, vareps, mll;
  double summisc = 0.0;

  gsl_vector_view Z = gsl_vector_view_array(sZ, m);
  gsl_vector * Z_cp = gsl_vector_alloc(m);
  gsl_matrix_view ZtZ = gsl_matrix_view_array(sZtZ, m, m);
  gsl_matrix * K = gsl_matrix_alloc(n, m);
  gsl_vector_view K_irow;
  gsl_matrix * r = gsl_matrix_alloc(n+1, m);
  gsl_vector_view r_row_t;
  gsl_vector_view r_row_tp1 = gsl_matrix_row(r, n);
  gsl_vector_set_zero(&r_row_tp1.vector);
  
  std::vector<gsl_matrix*> L(n);
  std::vector<gsl_matrix*> N(n+1);
  N.at(n) = gsl_matrix_calloc(m, m);
  
  int dim5cp = dim[5];
  
  KF_steady(dim, sy, sZ, sT, sH, 
    sR, V, sQ, sa0, sP0, 
    &mll, &v, &f, &invf, &vof, K, &L, tol, maxiter);

  convref = dim[5];
  dim[5] = dim5cp;
  
  if (convref == -1) {
    convref = n;    
  } else
    convref = ceil(convref * ksconvfactor[0]);
  nmconvref = n - convref;

  gsl_vector_view vaux;

  gsl_matrix * Mmm = gsl_matrix_alloc(m, m);

  gsl_matrix_view maux1, maux2;
  
  gsl_vector * var_eps = gsl_vector_alloc(n);
  vaux = gsl_vector_view_array(&f[0], n);
  gsl_vector_set_all(var_eps, msHsq);
  gsl_vector_div(var_eps, &vaux.vector);
  gsl_vector_add_constant(var_eps, *sH);

  for (i = n-1; i > -1; i--)
  {
    ip1 = i + 1;
    
    if (i != n-1)  //the case i=n-1 was initialized above
      r_row_tp1 = gsl_matrix_row(r, ip1);
    r_row_t = gsl_matrix_row(r, i);

    gsl_blas_dgemv(CblasTrans, 1.0, L.at(i), &r_row_tp1.vector, 
      0.0, &r_row_t.vector);
    gsl_vector_memcpy(Z_cp, &Z.vector);
    gsl_vector_scale(Z_cp, vof[i]);
    gsl_vector_add(&r_row_t.vector, Z_cp);

    N.at(i) = gsl_matrix_alloc(m, m);

    if (i < convref || i > nmconvref)
    {
      gsl_matrix_memcpy(N.at(i), &ZtZ.matrix);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, L.at(i), N.at(ip1), 0.0, Mmm);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Mmm, L.at(i), invf[i], N.at(i)); 
    } else {
      gsl_matrix_memcpy(N.at(i), N.at(ip1));
    }

    if (i < convref || i == nm1) {
      K_irow = gsl_matrix_row(K, i);
    }

    gsl_blas_ddot(&K_irow.vector, &r_row_tp1.vector, &epshat);

    epshat -= vof[i];
    epshat *= -*sH;

    if (i < convref || i > nmconvref)
    {
      maux1 = gsl_matrix_view_array(gsl_vector_ptr(&K_irow.vector, 0), 1, m);
      maux2 = gsl_matrix_view_array(gsl_vector_ptr(Z_cp, 0), 1, m);    
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &maux1.matrix, N.at(ip1),
        0.0, &maux2.matrix);

      vaux = gsl_vector_view_array(gsl_vector_ptr(var_eps, i), 1);
      gsl_blas_dgemv(CblasNoTrans, msHsq, &maux2.matrix, &K_irow.vector, 
        1.0, &vaux.vector);
      vareps = gsl_vector_get(&vaux.vector, 0);
    }

    summisc += pow(epshat, 2) + vareps; 

    gsl_matrix_free(L.at(i));    
    gsl_matrix_free(N.at(ip1));
  }

  res[0] = -0.5 * (n / sH[0] - summisc / sHsq);

  gsl_matrix_free(N.at(0));
  gsl_vector_free(Z_cp);
  gsl_vector_free(var_eps);
  gsl_matrix_free(r);
  gsl_matrix_free(K);
  gsl_matrix_free(Mmm);
}

void KFKSDS_steady_v2eta (int *dim, double *sy, double *sZ, double *sZtZ, 
  double *sT, double *sH, double *sR, 
  double *V, double *sQ, double *sa0, double *sP0, 
  double *tol, int *maxiter, double *ksconvfactor,
  double *res)
{
  int i, ip1, n = dim[0], m = dim[2], //ir = dim[3], 
    id = dim[6], convref, nmconvref;
  
  std::vector<double> v(n), f(n), invf(n), vof(n);
  //double v[n], f[n], invf[n], vof[n], Vsq = pow(V[0], 2),
  double Vsq = pow(V[0], 2), 
    etahat, vareta = 0.0, mll;
  double summisc = 0.0;

  gsl_vector_view Z = gsl_vector_view_array(sZ, m);
  gsl_vector * Z_cp = gsl_vector_alloc(m);
  gsl_matrix_view ZtZ = gsl_matrix_view_array(sZtZ, m, m);
  gsl_matrix * K = gsl_matrix_alloc(n, m);
  gsl_matrix * r = gsl_matrix_alloc(n+1, m);
  gsl_vector_view r_row_t;
  gsl_vector_view r_row_tp1 = gsl_matrix_row(r, n);
  gsl_vector_set_zero(&r_row_tp1.vector);
 
  std::vector<gsl_matrix*> L(n);
  std::vector<gsl_matrix*> N(n+1);
  N.at(n) = gsl_matrix_calloc(m, m);

  int dim5cp = dim[5];
  
  KF_steady(dim, sy, sZ, sT, sH, 
    sR, V, sQ, sa0, sP0, 
    &mll, &v, &f, &invf, &vof, K, &L, tol, maxiter);

  convref = dim[5];
  dim[5] = dim5cp;
  
  if (convref == -1) {
    convref = n;    
  } else
    convref = ceil(convref * ksconvfactor[0]);
  nmconvref = n - convref;

  gsl_matrix * Mmm = gsl_matrix_alloc(m, m);
  //gsl_matrix_view maux1, maux2;
  
  for (i = n-1; i > -1; i--)
  {
    ip1 = i + 1;
    
    if (i != n-1)  //the case i=n-1 was initialized above
      r_row_tp1 = gsl_matrix_row(r, ip1);
    r_row_t = gsl_matrix_row(r, i);

    gsl_blas_dgemv(CblasTrans, 1.0, L.at(i), &r_row_tp1.vector, 
      0.0, &r_row_t.vector);
    gsl_vector_memcpy(Z_cp, &Z.vector);
    gsl_vector_scale(Z_cp, vof[i]);
    gsl_vector_add(&r_row_t.vector, Z_cp);

    N.at(i) = gsl_matrix_alloc(m, m);
    if (i < convref || i > nmconvref)
    { 
      gsl_matrix_memcpy(N.at(i), &ZtZ.matrix);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, L.at(i), N.at(ip1), 0.0, Mmm);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Mmm, L.at(i), invf[i], N.at(i)); 
    } else {
      gsl_matrix_memcpy(N.at(i), N.at(ip1));
    }

    etahat = V[0] * gsl_matrix_get(r, ip1, id);

    if (i < convref || i > nmconvref)
    {
      vareta = V[0] - Vsq * gsl_matrix_get(N.at(ip1), id, id);
    }

    if (i != n-1)
    {
      summisc += pow(etahat, 2) + vareta;
    }

    gsl_matrix_free(L.at(i));    
    gsl_matrix_free(N.at(ip1));
  }

  res[0] = -0.5 * ((n-1) / V[0] - summisc / Vsq);
  
  gsl_matrix_free(N.at(0));
  gsl_vector_free(Z_cp);
  gsl_matrix_free(r);
  gsl_matrix_free(K);
  gsl_matrix_free(Mmm);
}
