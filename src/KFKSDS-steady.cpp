#include "KFKSDS-steady.h"

void KF_steady (int *dim, double *sy, double *sZ, double *sT, double *sH, 
  double *sR, double *sV, double *sQ, double *sa0, double *sP0, 
  double *mll, 
  std::vector<double> *v, std::vector<double> *f, 
  std::vector<double> *invf, std::vector<double> *vof,
  gsl_matrix *K, std::vector<gsl_matrix*> *L,
  double *tol, int *maxiter)
{
  int i, s, n = dim[0], p = dim[1], m = dim[2], r = dim[3], 
    conv = 0, counter = 0;
  double Kisum, Kim1sum;
  
  mll[0] = 0.0;

  // data and state space model matrices

  gsl_vector_view Z = gsl_vector_view_array(sZ, m);
  gsl_matrix_view T = gsl_matrix_view_array(sT, m, m);
  gsl_matrix_view Q = gsl_matrix_view_array(sQ, m, m);

  gsl_vector * a_pred = gsl_vector_alloc(m);
  //gsl_matrix * P_pred = gsl_matrix_alloc(m, m);
  
  // storage vectors and matrices

  gsl_vector * Vm = gsl_vector_alloc(m);
  gsl_vector * Vm_cp = gsl_vector_alloc(m);
  gsl_vector * Vm_cp2 = gsl_vector_alloc(m);
  gsl_matrix * Mmm = gsl_matrix_alloc(m, m);
  gsl_vector_view a0 = gsl_vector_view_array(sa0, m);
  gsl_vector * a_upd = gsl_vector_alloc(m);
  gsl_vector_memcpy(a_upd, &a0.vector);
  gsl_matrix_view P0 = gsl_matrix_view_array(sP0, m, m);
  gsl_matrix * P_upd = gsl_matrix_alloc(m, m);
  gsl_matrix_memcpy(P_upd, &P0.matrix);
  gsl_vector_view K_irow, K_im1row, Kri;
  gsl_matrix_view maux1, maux2;

  // filtering recursions

  for (i = 0; i < n; i++)
  {
    gsl_blas_dgemv(CblasNoTrans, 1.0, &T.matrix, a_upd, 0.0, a_pred);

    if (conv == 0) {
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &T.matrix, P_upd,
        0.0, Mmm);
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, Mmm, &T.matrix, 
        0.0, P_upd);
      gsl_matrix_add(P_upd, &Q.matrix);
    }

    gsl_blas_ddot(&Z.vector, a_pred, &v->at(i));
    v->at(i) = sy[i] - v->at(i);

    if (conv == 0) {    
      gsl_blas_dgemv(CblasNoTrans, 1.0, P_upd, &Z.vector, 0.0, Vm);
      gsl_blas_ddot(&Z.vector, Vm, &f->at(i));
      f->at(i) += *sH;
      invf->at(i) = 1.0 / f->at(i);
    } else {
      f->at(i) = f->at(i-1);
      invf->at(i) = invf->at(i-1);
    }

    gsl_vector_memcpy(Vm_cp, Vm);
    gsl_vector_memcpy(Vm_cp2, Vm);

    vof->at(i) = v->at(i) * invf->at(i);

    mll[0] += 0.5 * (log(M_2PI * f->at(i)) + (v->at(i) * vof->at(i)));
    
    if (conv == 0) {    
      maux1 = gsl_matrix_view_array(gsl_vector_ptr(Vm_cp, 0), m, 1);
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &maux1.matrix, 
        &maux1.matrix, 0.0, Mmm);
      gsl_matrix_scale(Mmm, invf->at(i));
      gsl_matrix_sub(P_upd, Mmm);
    }

    gsl_vector_memcpy(a_upd, a_pred);
    gsl_vector_scale(Vm_cp2, vof->at(i));
    gsl_vector_add(a_upd, Vm_cp2);

    K_irow = gsl_matrix_row(K, i);
    if (conv == 0) {    
      gsl_vector_scale(Vm_cp, invf->at(i));
      gsl_blas_dgemv(CblasNoTrans, 1.0, &T.matrix, Vm_cp, 0.0, &K_irow.vector);
    } else {
      K_im1row = gsl_matrix_row(K, i-1);
      gsl_vector_memcpy(&K_irow.vector, &K_im1row.vector);  
    }

    L[0].at(i) = gsl_matrix_alloc(m, m);
    if (conv == 0) {
      maux1 = gsl_matrix_view_array(gsl_vector_ptr(&K_irow.vector, 0), m, 1);
      maux2 = gsl_matrix_view_array(gsl_vector_ptr(&Z.vector, 0), 1, m);
      gsl_matrix_memcpy(L[0].at(i), &T.matrix);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, &maux1.matrix, 
        &maux2.matrix, 1.0, L[0].at(i));
    } else {
      gsl_matrix_memcpy(L[0].at(i), L[0].at(i-1));
    }

    // check if convergence to the steady state has been reached

    if (i > 0 & conv == 0)
    {
      if (i == 1)
      {
        Kim1sum = 0.0;
      }

      Kri = gsl_matrix_row(K, i);
      std::vector<double> vm((&Kri.vector)->data, (&Kri.vector)->data + m);
      Kisum = std::accumulate(vm.begin(), vm.end(), 0.0);

      if (fabs(f->at(i) - f->at(i-1)) < *tol)
      {
        counter += 1;
      }
      Kim1sum = Kisum;
      
      if (counter == *maxiter) {
        conv = 1;
        dim[5] = i;
      }
    }      
  }
  
  // deallocate memory
  
  gsl_vector_free(a_pred);
  gsl_vector_free(a_upd);
  gsl_matrix_free(P_upd);
  gsl_vector_free(Vm);
  gsl_vector_free(Vm_cp);
  gsl_vector_free(Vm_cp2);
  gsl_matrix_free(Mmm);
}

void KFKSDS_steady (int *dim, double *sy, double *sZ, double *sT, double *sH, 
  double *sR, double *sV, double *sQ, double *sa0, double *sP0,
  double *tol, int *maxiter, double *ksconvfactor,
  double *mll, double *epshat, double *vareps,
  double *etahat, double *vareta, 
  double *sumepsmisc, double *sumetamisc)
{
  int i, ip1, n = dim[0], m = dim[2], ir = dim[3], convref, nmconvref, nm1 = n-1;
  int irsod = ir * sizeof(double);

  //double v[n], f[n], invf[n], vof[n];
  std::vector<double> v(n), f(n), invf(n), vof(n);

  sumepsmisc[0] = 0.0;

  gsl_vector * sum_eta_misc = gsl_vector_calloc(ir);
  gsl_vector * etahat_sq = gsl_vector_alloc(ir);
  gsl_vector_view Z = gsl_vector_view_array(sZ, m);
  gsl_vector * Z_cp = gsl_vector_alloc(m);
  gsl_matrix * K = gsl_matrix_alloc(n, m);
  gsl_vector_view K_irow;
  gsl_matrix_view Q = gsl_matrix_view_array(sQ, m, m);
  gsl_matrix_view V = gsl_matrix_view_array(sV, ir, ir);  
  gsl_matrix_view R = gsl_matrix_view_array(sR, m, ir);

  gsl_matrix * r = gsl_matrix_alloc(n + 1, m);
  gsl_vector_view r_row_t;
  gsl_vector_view r_row_tp1 = gsl_matrix_row(r, n);
  gsl_vector_set_zero(&r_row_tp1.vector);

  std::vector<gsl_matrix*> L(n);
  std::vector<gsl_matrix*> N(n+1);
  N.at(n) = gsl_matrix_calloc(m, m);
  gsl_vector_view Ndiag;
  
  gsl_vector_view Qdiag = gsl_matrix_diagonal(&Q.matrix);
  gsl_vector * Qdiag_msq = gsl_vector_alloc(m);
  gsl_vector_memcpy(Qdiag_msq, &Qdiag.vector);
  gsl_vector_mul(Qdiag_msq, &Qdiag.vector);
  gsl_vector_scale(Qdiag_msq, -1.0);
  
  gsl_vector * sum_vareta = gsl_vector_calloc(m);

  KF_steady(dim, sy, sZ, sT, sH, 
    sR, sV, sQ, sa0, sP0, 
    mll, &v, &f, &invf, &vof, K, &L, tol, maxiter);

  convref = dim[5];
  if (convref == -1) {
    convref = n;    
  } else 
    convref = ceil(convref * ksconvfactor[0]);
  nmconvref = n - convref;

  gsl_vector_view vaux;

  gsl_matrix * Mmm = gsl_matrix_alloc(m, m);

  gsl_matrix * ZtZ = gsl_matrix_alloc(m, m);
  gsl_matrix_view maux1, maux2;
  maux1 = gsl_matrix_view_array(gsl_vector_ptr(&Z.vector, 0), m, 1);
  gsl_vector_memcpy(Z_cp, &Z.vector);
  maux2 = gsl_matrix_view_array(gsl_vector_ptr(Z_cp, 0), 1, m);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &maux1.matrix, 
    &maux2.matrix, 0.0, ZtZ);

  gsl_vector * var_eps = gsl_vector_alloc(n);

  double msHsq = -1.0 * pow(*sH, 2);
  vaux = gsl_vector_view_array(&f[0], n);
  gsl_vector_set_all(var_eps, msHsq);
  gsl_vector_div(var_eps, &vaux.vector);
  gsl_vector_add_constant(var_eps, *sH);

  gsl_matrix * eta_hat = gsl_matrix_alloc(n, ir);  
  gsl_matrix * Mrm = gsl_matrix_alloc(ir, m);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &V.matrix, &R.matrix, 0.0, Mrm);

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
      gsl_matrix_memcpy(N.at(i), ZtZ);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, L.at(i), N.at(ip1), 0.0, Mmm);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Mmm, L.at(i), invf[i], N.at(i)); 
    } else {
      gsl_matrix_memcpy(N.at(i), N.at(ip1));
    }
    
    if (dim[6] == 0 || dim[6] == 1)
    {

      if (i < convref || i == nm1) {
        K_irow = gsl_matrix_row(K, i);
      }

      gsl_blas_ddot(&K_irow.vector, &r_row_tp1.vector, &epshat[i]);

      epshat[i] -= vof[i];
      epshat[i] *= -*sH;

      if (i < convref || i > nmconvref)
      {
        maux1 = gsl_matrix_view_array(gsl_vector_ptr(&K_irow.vector, 0), 1, m);
        maux2 = gsl_matrix_view_array(gsl_vector_ptr(Z_cp, 0), 1, m);    
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &maux1.matrix, N.at(ip1),
          0.0, &maux2.matrix);

        vaux = gsl_vector_view_array(gsl_vector_ptr(var_eps, i), 1);
        gsl_blas_dgemv(CblasNoTrans, msHsq, &maux2.matrix, &K_irow.vector, 
          1.0, &vaux.vector);
        vareps[i] = gsl_vector_get(&vaux.vector, 0);
    } else {
        vareps[i] = vareps[ip1];
    }

    sumepsmisc[0] += epshat[i] * epshat[i] + vareps[i];
  }

  if (dim[6] == 0 || dim[6] == 2)
  {
    vaux = gsl_matrix_row(eta_hat, i);
    gsl_blas_dgemv(CblasNoTrans, 1.0, Mrm, &r_row_tp1.vector,
      0.0, &vaux.vector);

    memcpy(&etahat[i*ir], (&vaux.vector)->data, irsod);

    if (i != n-1)
    {
      gsl_vector_memcpy(etahat_sq, &vaux.vector);
      gsl_vector_mul(etahat_sq, etahat_sq);

      gsl_vector_add(sum_eta_misc, etahat_sq);
    }

    if (i != n-1)
    {
      if (i < convref || i > nmconvref)
      {
        Ndiag = gsl_matrix_diagonal(N.at(ip1));
        gsl_vector_memcpy(Z_cp, &Ndiag.vector);
        gsl_vector_mul(Z_cp, Qdiag_msq);
        gsl_vector_add(Z_cp, &Qdiag.vector);
        gsl_vector_set_zero(sum_vareta);
        gsl_vector_add(sum_vareta, Z_cp);
      }
        gsl_blas_dgemv(CblasTrans, 1.0, &R.matrix, sum_vareta, 1.0, sum_eta_misc);    
    }
  }

    gsl_matrix_free(L.at(i));    
    gsl_matrix_free(N.at(ip1));
  }

  gsl_matrix_free(N.at(0));

  if (dim[6] == 0 || dim[6] == 2)
  {
    memcpy(&sumetamisc[0], sum_eta_misc->data, irsod);
  }

  gsl_vector_free(Z_cp);
  gsl_vector_free(var_eps);
  gsl_vector_free(Qdiag_msq);
  gsl_vector_free(sum_vareta);
  gsl_vector_free(sum_eta_misc);
  gsl_vector_free(etahat_sq);
  gsl_matrix_free(eta_hat);  
  gsl_matrix_free(Mrm);
  gsl_matrix_free(r);
  gsl_matrix_free(K);
  gsl_matrix_free(ZtZ);
  gsl_matrix_free(Mmm);
}
