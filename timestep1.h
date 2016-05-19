// choose smallest and largest

odx1 = OARCL[2][1][k][j][i];
odx2 = OARCL[3][2][k][j][i];
if (odx1 < odx2) {
  ods = odx2;
  odl = odx1;
} else {
  ods = odx1;
  odl = odx2;
}

u = s[2][k][j][i];
rho = s[1][k][j][i];
// sound speed : used for idt[1]
if (wgam)
  cs2 = gam * (gam - 1.) * s[2][k][j][i] / s[1][k][j][i];
else
  cs2 = cs * cs;
cs = sqrt(cs2);			// need cs now


/* x-velocity plus or minus sound speed, as the limiting velocity */
vel1 = e2z_1(v[1][1][k], j, i) - vg[1];
if (vel1 >= 0)
{				// same as vel1=abs(vel1)+cs since sign 
				// not used below
vel1 = vel1 + cs;
} else {
  vel1 = vel1 - cs;
}

ftemp = invcour * vel1 * odx1;
idt2[2] = ftemp * ftemp;	// square takes care of sign

/* y-velocity plus or minus sound speed, as the limiting velocity */
vel2 = e2z_2(v[1][2][k], j, i) - vg[2];
if (vel2 >= 0)
{
vel2 = vel2 + cs;
} else {
  vel2 = vel2 - cs;
}

ftemp = invcour * vel2 * odx2;
idt2[3] = ftemp * ftemp;	// square takes care of sign


if (mag == 1) {
  /* alfven velocity */
  bxa = z2e_1(v[2][1][k], j, i);
  bya = z2e_2(v[2][2][k], j, i);
  valphen2 =
      (bxa * bxa + bya * bya +
       v[2][3][k][j][i] * v[2][3][k][j][i]) / rho;

  ftemp = invcour * ods;
  idt2[4] = ftemp * ftemp * valphen2;
} else {
  idt2[4] = SSMALL;
  valphen = 0.0;
  valphen2 = 0.0;
}

if (visc_art == 1) {
  /* linear viscosity */
#if(VISC_LINEAR)
  ftemp = invcour2 * nu_l * ods;
  idt2[5] = ftemp * ftemp * cs2;
#else
  idt2[5] = SSMALL;
#endif

  /* VNR viscosity: x-dir */
  dv = v[1][1][k][j][i + 1] - v[1][1][k][j][i];
  ftemp = invcour2 * nu_vnr * dv * odx1;
  idt2[6] = ftemp * ftemp;

  /* VNR viscosity: y-dir */
  dv = v[1][2][k][j + 1][i] - v[1][2][k][j][i];
  ftemp = invcour2 * nu_vnr * dv * odx2;
  idt2[7] = ftemp * ftemp;



} else {			// else if no art visc
  idt2[5] = SSMALL;
  idt2[6] = SSMALL;
  idt2[7] = SSMALL;
}

if (visc_real == 1) {

  ftemp = invcour2 * nu_real[k][j][i] * odx1 * odx1;
  idt2[8] = ftemp * ftemp;

#if(N2>1)
  ftemp = invcour2 * nu_real[k][j][i] * odx2 * odx2;
  idt2[9] = ftemp * ftemp;
#else
  idt2[9] = SSMALL;
#endif
} else {			// else if no real visc
  idt2[8] = SSMALL;
  idt2[9] = SSMALL;
}

if (res == 1) {
  /* resistivity */
  ftemp = invcour2 * res * ods;
  idt2[10] = ftemp * ftemp * cs2;
} else {
  idt2[10] = SSMALL;
}

